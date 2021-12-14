%function showWav(Y,Fs,varargin)
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
folderOut = getDirectory('fout');
buffer = 0;

Y = [];
Fs = 10000;
dura = 15;
showFig = 2;
time = datetime(2021,08,04,00,52,46);
%% Varagin
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                time = varargin{2};
            case 'file'
                file = vararing{2};
            case 'folderin'
                folderIn = varargin{2};
            case 'folderout'
                folderOut = varargin{2};
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

%% If no Y is specified try to read wav
if isempty(Y)
if ~exist('time'); error('Need to specify time.'); end
if ~exist('dura'); error('Need to specify duration.'); end
if ~exist('folderIn'); error('Need to specify wav location.'); end

% Read
[file, wavID] = getWavName(time, folderIn);
file = file{1};
fileTime = getFileTime(file);

[wav.s,wav.fs, time_s, audioInfo] = readBring([folderIn file], time,'duration',dura,'buffer',buffer);


convPW.SH =-194;      % Parameter for power convertion
convPW.G = 40;
convPW.D = 1;
Y =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;

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
    if isdatetime(time)
    itime = (second(time)-second(fileTime)) * Fs;    
    else
    itime = time * Fs;
    end
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


%%
for ii=1:nbSP
    [~,freq,time,tmp] = spectrogram(Y(:,i_ch(ii)),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    if ii==1; PSD2 = nan([nbSP size(tmp)]) ; end
    PSD2(ii,:,:) = 10*log10(tmp);
end
%%  ----------------------------------------------------------------
% -------------------- Figure 2 -----------------------------------
if any(showFig ==2)
    
% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 80];           % [dB] C limite pcolor
spgm.im.dur = dura;%'all';         % [s or 'all'] figure duration
%spgm.im.ovlp = 50;                 % [%] image window overlap
spgm.im.figvision = true ;          % [true false] visiblity of figure before saveas jpf file
% spectrogram window parameters
spgm.win.dur = 2048/10000;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 2;                  % [s] spectrogram zero padding
spgm.im.fmin = spgm.im.freqlims(1);spgm.im.fmax = spgm.im.freqlims(2);
spgm.fs = wav.fs; spgm.ns =  audioInfo.Ns; spgm.im.ns= spgm.ns; spgm.im.dur = audioInfo.dura; 
spgm = getSpgmWin(spgm);        % Get spectograme windows parameter     


% figure 1 : Show 10 chanel
sp.height=30; sp.width=30; sp.nbx=2; sp.nby=5;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.5; sp.spacey=1;
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');
%mesh2= mesh2';
nbSP = sp.nbx * sp.nby;
fontS = 12;

if ~exist('i_ch'); i_ch = 1:2:2*nbSP; end   

fig = figure(1);
movegui(fig,'northwest');
% --------------- Figure --------------------------
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [sp.width sp.height]);
set(fig, 'PaperPositionMode', 'manual');


i_ch = 1:2:2*nbSP;
for ii=1:nbSP
    [sy,sx] = find(mesh2==ii);
    % Get the spectogram of first chanel
    
    PSD = squeeze(PSD2(ii,:,:));
    ax(ii)=axes('position',sp.pos{sx,sy},'FontSize',fontS,'Box','on','Layer','top');
    pcolor(time,freq,PSD);
    shading flat;
    colormap(jet)
    %axis tight;
    ylim(spgm.im.freqlims);
    %set(gca,'YScale',spgm.im.freqscale);
    caxis([30 80])
    grid on
    
    
    % Text and label
    text(time(end) - 0.22 , spgm.im.freqlims(2) * 0.93,['ch.' num2str(i_ch(ii))], 'backgroundcolor', [1 1 1],'margin',2)
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2 + 5,['ch.'])
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2-5,[num2str(mainCh(ii))])
    if sx == 1
        yl = ylabel('f (Hz)','FontSize',fontS);
    else
        ax(ii).YTickLabel = [];
    end
    if sy == sp.nby
        xl = xlabel('t (s)','FontSize',fontS) ;
    else
        ax(ii).XTickLabel = [];
    end
    
end

% Color bar
cbPos = [sp.pos{sx,sy}(1)+sp.pos{sx,sy}(3)+0.01 sp.pos{sx,sy}(2) 0.02 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
%ax(end).Position = sp.pos{end};
ylabel(cb,'PSD (dB re. 1µPa2/Hz)')

% General text and label
sgtitle([arrID ' | ' datestr(ptime)],'Interpreter','none');

print('-dpng','-r150',[folderOut 'showWavChanels' moreInfo  '.png'])

end
