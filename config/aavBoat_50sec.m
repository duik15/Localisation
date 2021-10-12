    % Config parameter for locateBring.m
    addpath('../')
    
    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'AAV';
    outName = 'aavCircleBoat_Ns14_f150-200hz';
    folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];  
    pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircle_7sec_Ns16_f150-200hz/'];
    openData = true;
    
    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
    
    % Loading files and time
    % Time must be in datetime format. The file to load will be automatically find
    ptimeB = getPingInfo(arrID);
    [ptimeA plocA] = loadGPXTrack([],arrID);
    % Make a new time from start to end of waht I have
    ifirst = find(  0<(ptimeA - min(ptimeB)), 1 , 'first' ); 
    iend = find(  0<(ptimeA - max(ptimeB)), 1 , 'first' ); 
    ptime = ptimeA(ifirst:iend);
    ploc = plocA(ifirst:iend,:);
    
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
    arrLoc = getArrLoc(arrID);
    angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));
    
    
  
    % If wanted to open alredy run script
    if openData == true
        disp('Opening already run data')
        p = getRunData(pingFolder)
        bName = dir([folderOut 'data*']);
        load ([folderOut bName.name])
        
        % Ploting the global figure
        showFig =10;
        showFigAngleTime;
    end