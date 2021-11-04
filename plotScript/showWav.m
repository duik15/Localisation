function showWav(Y,Fs,varargin)
% showWav(Y,Fs)
% showWav(Y,Fs,'time', 120, 'duration',20,'chanel',[1 2 3],'tPos',1,'printFig', true)
disp('< -- Automatic wav plot -- >')

% Default Parameters
showFig = [1];   % Plot ID, can be an array for mutiple plot
moreInfo = '';  % Add info when saving XX_moreInfo.png
%Fs = 10000;    % Sample rate of the file
chanel = 1;     % The chanel to plot
tPos = 1 ;      % Time windows position when specifing the time. 1= start, 2= middle 3= end
%dura = 20;  % Lenght of time axis in second
resolution = '-r150';   % Figure resolution
printFig = false;

%% Varagin
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                time = varargin{2};
            case 'showfig'
                showFig = varargin{2}; 
            case 'chanel'
                chanel = varargin{2};
            case 'moreinfo'
                moreInfo = varargin{2}; 
            case 'duration'
                dura = varargin{2};
            case 'climit'
                CLimit = varargin{2};
            case 'fs'
                Fs = varargin{2};
            case 'tpos'
                tPos = varargin{2};
            case 'resolution'
                resolution = varargin{2};
            case 'printFig'
                printFig = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end
%% Automatic parameter
sY = size(Y);
if any(chanel > sY(2))
       disp('Chanel requested is too high.') 
end

% Taking only one chanel
Y1 = Y(:,chanel(1));   
   
% Setting the leght of the time axis
if ~exist('dura')
   dura  = length(Y1)/Fs; 
end

 % Loading the time iteration
if ~exist('time')    
    if tPos ==1;        itime= 1 ; 
    elseif tPos == 2;   itime= length(Y1)/2 ; 
    else                itime= length(Y1);
    end
else
    itime = time * Fs;
end

% Start and end point
if tPos ==1
    istart = itime;
    iend = itime + Fs*dura-1;
elseif tPos ==2
    istart = itime - Fs*dura/2;
    iend = itime + Fs*dura/2 -1;
else
    istart = itime;
    iend = itime - Fs*dura -1 ;
end

%% Ploting figure
if any(showFig==1)

figure(1)
plot(Y1(istart:iend))
ax = gca;

ax.XTickLabel = ax.XTick /Fs;

ylabel('Amplitude')
xlabel('Time (s)')

if printFig == true
    print('-dpng','-r150',['showWav' moreInfo  '.png']);
end
end