clear all
close all

% Config parameter for locateBring.m
addpath(genpath('../'));

% Need to run data or just open it?
openData = false ;

%repertoire des fonctions m_map
%repmap = 'T:\Utilitaires\Matlab_tools\m_map2020\m_map\';
%addpath(genpath(repmap))
% Path information : folderIn = wav folder / folderOut = figure output folder
% arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
arrID = 'CLD';
outName = 'ctd5_Ns14_f150-200hz';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];
folderIn = ['/Volumes/KevDD/BRING/wav/CLD/'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];
%boatFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/mlbBoatTrack_Ns14_f150-200hz/'];

% Loadign data
ptime= datetime(2021,07,14,01,58,15);
ploc=[49.2237, -64.6179];
[aloc aoff arrOri]= getArrLoc(arrID);
[angleR distR]= getRealAngle(arrID, ploc(:,1),ploc(:,2));
nbP = length(ploc);

% Figure parameters
showFig = [1 2 3]       % Figure number to print
saveData = true;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^14;              % Total number of sample
buffer = 0.4;             % Time in second

% Spectro parameter
% Frequence min and max
fmin_int = 100;
fmax_int = 150;

% Get the real angle
arrLoc = getArrLoc(arrID);
angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));


% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    b = getRunData(boatFolder);
else
    locateBring;
end


if ~isfolder(folderOut ); disp(['Creating output folder: ' outName]); mkdir(folderOut); end
%% Create fake data
%[lon,lat, phi] = m_fdist(aloc(2),aloc(1), ((ii-1)* 360 / nbH ) + theta , r);

%for i=1:nbP
%    [distap(i,1),aap(i,1),apa2(i,1)] = m_idist(aloc(2),aloc(1),ploc(i,2),ploc(i,1));
%    apa(i,1) = apa2(i) + 10 *(2* rand() - 1);
%end

%% Figure 1
close all
% Zone de la carte
winSize = 0.015;
%winSize = 0.15;
lon_min = -aloc(2) - winSize ;
lon_max = -aloc(2) + winSize;
lat_min =  aloc(1) -  winSize/2;
lat_max =  aloc(1) + winSize/2;

% Les isobathes à contourer
isobath = [0:20:200];

% Dimensions de la figure sur papier (en cm)
paperwidth  = 12;
paperheight = 12;


figure(1)
clf
set(groot, 'defaultAxesTickLabelInterpreter','tex');
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('box','fancy','fontsize',8,'linestyle','none','fontname','times');
%m_grid('box','fancy','fontsize',8,'linestyle','none','fontname','times','xticklabel',{'49°W'});

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco_colormap);

hold on
m_contourf(LON,LAT,Z,isobath,'linestyle','none');
% Information graphique
c = colorbar('location','eastoutside');
%set(c,'Position',[0.2 0.2 0.4 0.02],'FontSize', 8, 'FontName','times','xAxisLocation','top');
cpos = c.Position;
set(c, 'FontSize', 8, 'FontName','times')
ylabel(c,'Depth(m)')
%m_text(-69.7,46.9,'Depth (m)','FontSize', 8, 'FontName','times','HorizontalAlignment','center');
xlabel('Longitude', 'FontSize', 10, 'FontName','times')
ylabel('Latitude', 'FontSize', 10, 'FontName','times')
set (gca,'TicklabelInterpreter','tex')

m_gshhs_l('patch',[.7 .7 .7]); % coastlines



% Hydrophone location
hc = m_plot(-aloc(2),aloc(1),'+','MarkerSize',14);

% Circular ping
for i=1:length(ploc)
    hp1=   m_plot(ploc(i,2),ploc(i,1),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','green','MarkerSize',4) ;
    %m_text(-ploc(i,2) + 0.0001 ,ploc(i,1) + 0.0001 ,num2str(i),'FontSize', 14);
    %m_text(-p.ploc(i,2) + 0.0001 ,p.ploc(i,1) + 0.0001 ,num2str(p.angleM()) ,'FontSize', 14);
end


% Now lets find our location
apa = p.angleM + 180;
% Circular ping
for i=1:length(ploc)
    % Calculate the orietation line
    [lonin,latin, phi] = m_fdist(ploc(i,2),ploc(i,1), 360- apa(i), 700);
    [loninc(i),latinc(i), phi] = m_fdist(ploc(i,2),ploc(i,1), 360 - ( apa(i)+180), 500);
    %hl = m_plot(-[ploc(i,2) lonin], [ploc(i,1) latin],'-.k');
    hl = m_plot([ploc(i,2) lonin], [ploc(i,1) latin],'-.k');
    %m_text(-ploc(i,2) + 0.0001 ,ploc(i,1) + 0.0001 ,num2str(i),'FontSize', 14);
end


if exist('b')
    %Now lets find our location
    apab = b.angleM + 180;
    % Circular ping
    for i=1:length(b.ploc)
        % Calculate the orietation line
        [loninb,latinb, phib] = m_fdist(-b.ploc(i,2),b.ploc(i,1), 360- apab(i), 700);
        [lonincb(i),latincb(i), phibb] = m_fdist(-b.ploc(i,2),b.ploc(i,1), 360 - ( apab(i)+180), 500);
    end
end

%hh = m_plot(-mean(ploc(:,2)),mean(ploc(:,1)),'x','MarkerSize',14);
hc2 = m_plot(mean(loninc),mean(latinc),'xr','MarkerSize',14);
%hb = m_plot(mean(lonincb),mean(latincb),'xg','MarkerSize',14);


% Print legend location
hct = m_plot(-aloc(2)- 0.013,aloc(1) + 0.006,'+','MarkerSize',14,'Color',hc.Color);
m_text(-aloc(2)- 0.012,aloc(1) + 0.006,[num2str(aloc(1)) ' N, ' num2str(aloc(2)) ' W'],'FontSize', 14);
hct2 = m_plot(-aloc(2)- 0.013,aloc(1) + 0.005,'x','MarkerSize',14,'Color',hc2.Color);
m_text(-aloc(2)- 0.012,aloc(1) + 0.005,[num2str(mean(latinc)) ' N, ' num2str(mean(loninc)) ' W'],'FontSize', 14);


% Legend
leg = m_legend([hc,hp1,hc2], 'Circular Array','500m Ring','Retreived center')
%leg.Position = [49,-64,10,4]
leg.Position = [49,5,25,3]


%%
print('-dpng', '-r300',[folderOut outName 'crossLineCircle.png']);

