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
    ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
    
    % Figure parameters
    showFig = [1 2 3 4 7]       % Figure number to print
    saveData = true;   % Save result to .mat
    printFig = true;    % Saving figure to a folder
    nbPk = 4 ;          % Nomber of side lobe to keep
    
    % Reading parameters
    Ns = 2^16;              % Total number of sample
    buffer = 0.4;             % Time in second
    
    
    % Spectro parameter
    % Frequence min and max
    fmin_int = 150;
    fmax_int = 200;
    
     % Get the real angle
    arrLoc = getArrLoc(arrID);
    angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));