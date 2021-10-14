% This script show the location of ping emission

%repertoire des fonctions m_map
%repmap = 'T:\Utilitaires\Matlab_tools\m_map2020\m_map\';
%addpath(genpath(repmap))

% Parameter
locID = 'av';

% Loadign data
[ploc, ptime] = getPingLoc(['c' locID]);
aloc = getArrLoc(locID);


% Read bathymetry
load ~/Documents/MATLAB/Oceanographie/m_map/mmap_ex/GSL_bathy_500m.mat;
load ~/Documents/MATLAB/Oceanographie/m_map/mmap_ex/gebco_colormap.dat;


%% Figure 1

% Zone de la carte
winSize = 0.06;
%winSize = 0.15;
lon_min = -aloc(2) - winSize ;
lon_max = -aloc(2) + winSize;
lat_min =  aloc(1) -  winSize/2;
lat_max =  aloc(1) + winSize/2;

% Les isobathes à contourer
isobath = [0:10:450];

% Dimensions de la figure sur papier (en cm) 
paperwidth  = 12;
paperheight = 12;


figure(1)
clf

m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('box','fancy','fontsize',8,'linestyle','none','fontname','times');

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
set (gca,'TicklabelInterpreter','None')

m_gshhs_f('patch',[.7 .7 .7]); % coastlines
  
% Printing location

% Hydrophone location
hh = m_plot(-aloc(2),aloc(1),'p','MarkerEdgeColor','k',...
    'MarkerFaceColor','red','MarkerSize',14);
%m_text(-aloc(2) + 0.01 ,aloc(1) +                                          0.01,nameStat{1},'FontSize', 14);


% Circular ping
for i=1:length(ploc)
hp1=   m_plot(-ploc(i,2),ploc(i,1),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','green','MarkerSize',4) ;
%m_text(-ploc(2) + 0.01 ,ploc(1) + 0.01,,'FontSize', 14);
end

% Legend
leg = m_legend([hh,hp1], 'Circular Array','500m Ring')
leg.Position = [50,64,10,4]

print('-dpng', '-r300',['results/AV_real_circularPing.png']);    

