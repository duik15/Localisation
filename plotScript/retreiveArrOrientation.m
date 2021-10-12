% This script show the location of ping emission
clear all
close all

%repertoire des fonctions m_map
%repmap = 'T:\Utilitaires\Matlab_tools\m_map2020\m_map\';
%addpath(genpath(repmap))
outName = getOutName('retreive')
moreInfo = 'av_reverse';

% Parameter
locID = 'av';

% Loadign data
[ploc, ptime] = getPingLoc(['c' locID]);
aloc = getArrLoc(locID);
nbP = length(ploc);

% Read bathymetry
load ~/Documents/MATLAB/Oceanographie/m_map/mmap_ex/GSL_bathy_500m.mat;
load ~/Documents/MATLAB/Oceanographie/m_map/mmap_ex/gebco_colormap.dat;

% Load data
load '/Users/Administrator/Documents/MPO/BRing/Data/results/AV/realOrder_OffSet23_voc3_Ns2e14_AvCircle/locateBRING_Data_realOrder_OffSet23_voc3_Ns2e14_AvCircle.mat'
apa = aMax -90;

if ~isfolder(outName); disp(['Creating output folder: ' outName]); mkdir(outName); end
%% Create fake data
%[lon,lat, phi] = m_fdist(aloc(2),aloc(1), ((ii-1)* 360 / nbH ) + theta , r);

%for i=1:nbP
%    [distap(i,1),aap(i,1),apa2(i,1)] = m_idist(aloc(2),aloc(1),ploc(i,2),ploc(i,1));
%    apa(i,1) = apa2(i) + 10 *(2* rand() - 1);
%end

%% Figure 1
close all
% Zone de la carte
winSize = 0.02;
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
hh = m_plot(-aloc(2),aloc(1),'+','MarkerSize',14);
%m_text(-aloc(2) + 0.01 ,aloc(1) + 0.01,nameStat{1},'FontSize', 14);


% Circular ping
for i=1:length(ploc)
hp1=   m_plot(-ploc(i,2),ploc(i,1),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','green','MarkerSize',4) ;
%m_text(-ploc(i,2) + 0.0001 ,ploc(i,1) + 0.0001 ,num2str(i),'FontSize', 14);
m_text(-ploc(i,2) + 0.0001 ,ploc(i,1) + 0.0001 ,num2str(apa(i)+180) ,'FontSize', 14);
end

% Now lets find our location
apa2 = apa -180;
% Circular ping
for i=1:length(ploc)
%i=1;

% Calculate the orietation line
[lonin,latin, phi] = m_fdist(ploc(i,2),ploc(i,1), apa(i) +180, 700);
[loninc(i),latinc(i), phi] = m_fdist(ploc(i,2),ploc(i,1), apa(i)+180, 500);
%hl = m_plot(-[ploc(i,2) lonin], [ploc(i,1) latin],'-.k');
hl = m_plot(-[ploc(i,2) lonin], [ploc(i,1) latin],'-.k');
%m_text(-ploc(i,2) + 0.0001 ,ploc(i,1) + 0.0001 ,num2str(i),'FontSize', 14);
end

%hh = m_plot(-mean(ploc(:,2)),mean(ploc(:,1)),'x','MarkerSize',14);
hc = m_plot(-mean(loninc),mean(latinc),'xr','MarkerSize',14);

% Legend
leg = m_legend([hh,hp1,hc], 'Circular Array','500m Ring','Retreived center')
%leg.Position = [49,-64,10,4]
leg.Position = [49,5,25,3]
%%
print('-dpng', '-r300',[outName 'crossLineCircle' moreInfo '.png']);    

