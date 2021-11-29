    % Config parameter for locateBring.m
    addpath('../')

    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'AAV';
    outName = 'charlotte_aavCircle_7sec_Ns16_f150-200hz';
    folderIn = ['E:\Bring_Dep_1\' arrID '\']; % Local Mac folder
    folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];    
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
    
    
    
    % If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    %unpackStruct(p);  % Unpack structure to get same workspace as locateBring
else
    locateBRing; % The main loop calculation are locate in this script
    showGlobalFig;
    % Add some figure script located in ../plotScript  
    %showAOACircle; 
end