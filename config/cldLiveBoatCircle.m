    % Config parameter for locateBring.m
    addpath(genpath('../'));
    
    % Need to run data or just open it?
    openData = true;
    
    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'CLD';
    outName = 'CLDLiveBoat_onlyTrackTime_Ns14_f150-200hz';
    folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];  
    pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircle_7sec_Ns16_f150-200hz/'];
    boatFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircleBoat_Ns14_f150-200hz/'];
  
    
    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
    
    % Loading files and time
    % Time must be in datetime format. The file to load will be automatically find
    [ptimeT plocT] = loadGPXTrack([],arrID);
    ptimeP = getPingInfo(arrID);
    [indexT val ] = getInbetween(ptimeT, [min(ptimeP) max(ptimeP)]);

    % Make a new time from start to end of waht I have
    ptimeB = ptimeT(indexT);
    plocB = plocT(indexT,:);
    
    % Now opening file with a dt specifiy
    dt = seconds(4);
    ptime = min(ptimeB):dt:max(ptimeB);
    
    
    % Figure parameters
    showFig = []       % Figure number to print
    saveData = true;   % Save result to .mat
    printFig = false;    % Saving figure to a folder
    nbPk = 4 ;          % Nomber of side lobe to keep
    
    % Reading parameters
    Ns = 2^14;              % Total number of sample
    buffer = 0.4;             % Time in second
    
    
    % Spectro parameter
    % Frequence min and max
    fmin_int = 150;
    fmax_int = 200;
    
    % Get the real angle
    %arrLoc = getArrLoc(arrID);
    %angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));
    
    
  
    % If wanted to open alredy run script
    if openData == true
        disp('Opening already run data')
        p = getRunData(pingFolder);
        b = getRunData(boatFolder);
        bName = dir([folderOut 'data*']);
        load ([folderOut bName.name])
        
        % Ploting the global figure
        showFig =[10 11];
        showFigAngleTime;
    end