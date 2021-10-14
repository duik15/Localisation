    % Tracking the boat location at each gpx time. 
    clear all
    close all
    % Config parameter for locateBring.m
    addpath(genpath('../'));
    
    % Need to run data or just open it?
    openData = true;
    
    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'CLD';
    outName = 'CLDBoatTrack_Ns14_f150-200hz';
    folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];  
    pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/CLD/cldCircle_7sec_Ns16_f150-200hz/'];
    boatFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/CLD/CLDBoatTrack_Ns14_f150-200hz/'];
  
    
    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
    
    % Loading files and time
    % Time must be in datetime format. The file to load will be automatically find
    [ptimeT plocT] = loadGPXTrack([],arrID);
    ptimeP = getPingInfo(arrID);
    [indexT val ] = getInbetween(ptimeT, [min(ptimeP) max(ptimeP)]);

    % Make a new time from start to end of waht I have
    ptime = ptimeT(indexT);
    ploc = plocT(indexT,:);
    
    
    
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
        p = getRunData(pingFolder);
        b = getRunData(boatFolder);
        unpackStrcut(b);
        %bName = dir([folderOut 'data*']);
        %load ([folderOut bName.name])
        showFig = 10;
        showFigAngleTime;
        
    else
        locateBring;
    end
    
    %% ---------------------- Figure -----------------------------
    showpColorAngleTime;
    
    
    % pcolormap of the matrice energie
figure(2)

ax=gca;
pcolor(ptime, 1:360, 10*log10(matEnergie'))
shading flat
cb = colorbar('location','eastoutside');
ylabel(cb, 'Energy (dB)');
colormap jet

if exist('p')
hold on
hpm = plot(p.ptime,p.angleM,'p','color','r','markersize',15,'markerfacecolor','r');
end


% Axes paremeter
ylabel('Azimut $^{\circ}$','interpreter','latex')
xlabel('Time')

leg = legend([hpm hbm],'Ping azimut','Real boat azimut','Boat azimut');
%leg.Box ='off';

print('-dpng','-r150',[folderOut 'pcolorAngleTime_' outName '.png' ])