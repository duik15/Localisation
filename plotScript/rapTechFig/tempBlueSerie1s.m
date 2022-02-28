% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette
clear all
close all
%addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
%arrID = 'AAV'; 
arrID = 'CLD';


% Selected time
%ptime= datetime(2021,07,16, 06, 19 , 03 );
%ptime= datetime(2021,07,13, 23, 45, 10);

dList = [datetime(2021,07,13, 23, 45, 13) datetime(2021,07,13, 23, 45, 16);
               datetime(2021,07,13, 23, 45, 31) datetime(2021,07,13, 23, 45, 35);
               datetime(2021,07,13, 23, 45, 45) datetime(2021,07,13, 23, 45, 49);
               datetime(2021,07,13, 23, 45, 59) datetime(2021,07,13, 23, 41, 04);
               datetime(2021,07,13, 23, 46, 17) datetime(2021,07,13, 23, 46, 22);
               datetime(2021,07,13, 23, 46, 32) datetime(2021,07,13, 23, 46, 37);
               datetime(2021,07,13, 23, 46, 52) datetime(2021,07,13, 23, 46, 57);
               datetime(2021,07,13, 23, 47, 05) datetime(2021,07,13, 23, 47, 09);
               datetime(2021,07,13, 23, 47, 32) datetime(2021,07,13, 23, 47, 37);
               datetime(2021,07,13, 23, 48, 02) datetime(2021,07,13, 23, 48, 06);
               datetime(2021,07,13, 23, 48, 19) datetime(2021,07,13, 23, 48, 23);
               datetime(2021,07,13, 23, 48, 34) datetime(2021,07,13, 23, 48, 40);
               datetime(2021,07,13, 23, 48, 50) datetime(2021,07,13, 23, 48, 55);
               datetime(2021,07,13, 23, 49, 29) datetime(2021,07,13, 23, 49, 34);
               datetime(2021,07,13, 23, 49, 52) datetime(2021,07,13, 23, 49, 56)];

% Path information | You can also add and use fout = getDirectory('fout')
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderIn = ['/Volumes/BringDD4/Bring_Dep_1/' arrID '/']; % Local Mac folder
%folderIn = ['/Volumes/BringDD3/Bring_Dep_1/' arrID '/']; % Local Mac folder
folderIn = ['/Volumes/BringDD4/Bring_Dep_1/' arrID '/']; % Local Mac folder

%outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
%outName = [arrID '_' wavi.wavID '_mom' ];%'MLB_1493_20210804T005254';
outName = 'blueSerie';
%folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];
%%
% Figure parameters
showFig = [1 2 4 5]%2 3 4 5 ];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 2;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = imDur/2 ;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% Creating the output folder
%if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [20 100];       % [Hz] frequency scale boundary limits
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
close all
%for ii=1:size(dList,2)
    %ii=4 %AAV 
    ii=7 %CLD
    if strcmp(arrID,'AAV')
        indArr = 1;
    else
        indArr = 2;
    end
    
    % Get the file info and folder ready
    ptime = dList(ii,indArr)
    [~, wavi] = getWavName(ptime, folderIn);
    outName2 = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' outName2 '/'];
    if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
    
    disp(['Executing beamforming for one time at ' datestr(ptime)])
    mainBring;

    % Move folder name to add azimut
    if imDur <= 10
    movefile(folderOut,[folderOut(1:end-1) '_' num2str(angleM) 'd'])
    else
    movefile(folderOut,[folderOut(1:end-1) '_' num2str(imDur) 's'])
    end
    
    
%end
%% 



%% Add some more specified line or figures related to you run
% Enjoy!

disp('Done !')
