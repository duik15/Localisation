% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette
clear all
close all
%addpath(genpath('../'));

% Exploration vector
load(['/Users/Administrator/Documents/MPO/BRing/Data/wav/' 'DETECT_NARW_PRC_PR_V1.mat'])
nbF = numel(dtime);

% Main loop parameter
arrID = 'PRC';
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderIn = ['/Volumes/BringDD4/Bring_Dep_2/' arrID '/']; % Local Mac folder
%pathOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/explo/'];
folderIn = ['\\169.254.47.215\usbshare1-2\Bring_Dep_2\' arrID '\']
pathOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\explo\'];



% Figure parameters
showFig = [1 2 3 4 5 6];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 3;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = 1;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

%% Spectrogram parameter
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

%%
for ii=2:2%:nbF
disp(['Exploring file' getWavName(dtime(ii),folderIn)]) 
    
% Selected time
ptime = dtime(ii);
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
folderOut = [pathOut outName '/' ];

% Creating the output folder
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Executing 
mainBring;    


end % Of big loop
disp('Done !')
