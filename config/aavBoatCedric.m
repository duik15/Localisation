    clear all
    close all
    % Config parameter for locateBring.m
    addpath('../')
    
    % Need to run data or just open it?
    openData = false;
    
    % Path information : folderIn = wav folder / folderOut = figure output folder
    % arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
    arrID = 'AAV';
    outName = 'aavBoatCedric_onlyTrackTime_Ns14_f150-200hz';
    %folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
    folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];  
    %pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircle_7sec_Ns16_f150-200hz/'];
    %boatFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircleBoat_Ns14_f150-200hz/'];
    folderIn = ['/Volumes/KevDD/BRING/wav/' arrID '/']; % Local Mac folder
    
    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
    
    % Loading files and time
    % Time must be in datetime format. The file to load will be automatically find
    dt = seconds(4);
    t1 = datetime(2021,07,15,00,30,00);
    t2 = datetime(2021,07,15,00,45,00);
    ptime = t1:dt:t2;
    
    
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
    else
        locateBring2;
    end
    
    
    %%  ---------------------- figure ---------------------------
    
    % pcolormap of the matrice energie
figure(2)

ax=gca;
pcolor(ptime, 1:360, 10*log10(matEnergie'))
shading flat
cb = colorbar('location','eastoutside');
ylabel(cb, 'Energy (dB)');
colormap jet


hold on
hm = plot(ptime,angleM,'o','color','k','markersize',5);



% Axes paremeter
ylabel('Azimut $^{\circ}$','interpreter','latex')
xlabel('Time')

leg = legend([hm],'First  Lobe','Location','northwest');
%leg.Box ='off';

print('-dpng','-r150',[folderOut 'pcolorAngleTime_' outName '.png' ])
%{
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
if exist('b')
hold on
hbr = plot(b.ptime,b.angleR,'o','color','k','markersize',5,'markerfacecolor','k');
hbm = plot(b.ptime,b.angleM,'o','color','b','markersize',5,'markerfacecolor','b');
end

% Axes paremeter
ylabel('Azimut $^{\circ}$','interpreter','latex')
xlabel('Time')

leg = legend([hpm hbr hbm],'Ping azimut','Real boat azimut','Boat azimut');
%leg.Box ='off';

print('-dpng','-r150',[folderOut 'pcolorAngleTime_' outName '.png' ])
    %}