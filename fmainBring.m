function [angleM res] = fmainBring(arrID, ptime, folderIn ,varargin)

% Figure parameters
showFig = 0;       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 3;              % Total number of sample
buffer = imDur/2;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);

% Folder to print result
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
folderOut = [ getDirectory('fout')  arrID '/' outName '/' ];


% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 80];           % [dB] C limite pcolor
spgm.im.dur = imDur;%'all';         % [s or 'all'] figure duration
%spgm.im.ovlp = 50;                 % [%] image window overlap
spgm.im.figvision = true ;          % [true false] visiblity of figure before saveas jpf file
% spectrogram window parameters
spgm.win.dur = 2048/10000;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 2;                  % [s] spectrogram zero padding
spgm.im.fmin = spgm.im.freqlims(1);spgm.im.fmax = spgm.im.freqlims(2);

%% Input parameter reading

while ~isempty(varargin) 
    switch lower(varargin)
        case 'showfig'
            showFig = varagin{2};
        case 'imdur'
            imDur = varagin{2};
        case 'buffer'
            buffer = varagin{2};
        case 'spgm'
            spgm = varagin{2};
        case 'folderout'
            folderOut = varagin{2};           
        otherwise
            error(['Can''t understang input parameter ' varargin{1} ])
    end
    varargin(1:2)=[];
end


%%
% Creating the output folder
if any(showFig) && ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Execute the mainBring routine
mainBring;

% Save result to structure
res.angleM = angleM;
res.angleA = angleA;
res.dtime = ptime;
res.file = fileList{1};
res.wavID = wavID;
res.energie_dB = Energie_dB;
end