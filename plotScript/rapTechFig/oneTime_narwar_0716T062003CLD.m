% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette
clear all
close all
%addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
%arrID = 'CLD';
arrID = 'AAV';

% Selected time
%ptime= datetime(2021,07,16, 06, 19 , 03 );
ptime= datetime(2021,07,16, 06, 18, 56);

% Path information | You can also add and use fout = getDirectory('fout')
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderIn = ['/Volumes/BringDD4/Bring_Dep_1/' arrID '/']; % Local Mac folder
folderIn = ['/Volumes/BringDD3/Bring_Dep_1/' arrID '/']; % Local Mac folder
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
%outName = [arrID '_' wavi.wavID '_mom' ];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];

% Figure parameters
showFig = [1 2 3 4 5 6]%2 3 4 5 ];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  354;        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 2;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = imDur/2 - 0.2;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% Creating the output folder
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 60];           % [dB] C limite pcolor
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

% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    %unpackStruct(p);  % Unpack structure to get same workspace as locateBring
else
    % For figure related to locateBring see plotScript/showFigBring.m
    disp(['Executing beamforming for one time at ' datestr(ptime)])
    mainBring;
    %locateBRing;
    
    % Add some figure script located in ../plotScript
    % Show global figure
    %showGlobalFig;
    %showAOACircle; 
    
end
%% 

% Move folder name to add azimut
if imDur <= 4
movefile(folderOut,[folderOut(1:end-1) '_' num2str(angleM) 'd'])
else
    movefile(folderOut,[folderOut(1:end-1) '_' num2str(imDur) 's'])
end

%% Add some more specified line or figures related to you run
% Enjoy!

disp('Done !')
