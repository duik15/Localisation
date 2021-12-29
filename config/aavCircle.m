    % Config parameter for locateBring.m
    addpath('../')

    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'AAV';
    outName = 'aavCircle_7sec_Ns16_f150-200hz';
    folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];    
    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
    
    % Loading files and time
    % Time must be in datetime format. The file to load will be automatically find
    [ptime, ploc]  = getPingInfo(arrID); % Load ping information
    ptime = ptime + seconds(5);          % Add an offset to be center the 5 upcalls
    
    % Figure parameters
    showFig = [1 2 3 4 7]       % Figure number to print
    saveData = true;   % Save result to .mat
    printFig = true;    % Saving figure to a folder
    nbPk = 4 ;          % Nomber of side lobe to keep
    
    % Reading parameters
    %Ns = 2^16;              % Total number of sample
    buffer = 0.4;             % Time in second
    imDur = 4;
    
    % Spectro parameter
    % Frequence min and max
    fmin_int = 150;
    fmax_int = 200;
    
     % Get the real angle
    arrLoc = getArrInfo(arrID);
    angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));
    
    
    %% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [75 225];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 100];           % [dB] C limite pcolor
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

mainBring;