% This script show the 5 panel diagnostic for a double localisation
% Blue CLD_0950: 19P, 1:56L, 2:28P, 2:44L, 4:41L
% Blue AAV_0950: 16P, 1:56L*, 2:26P, 4:39P
clear all
close all

moreInfo = 'blueHonguedo';

arrID = {'CLD','AAV'};  
dtime = [datetime(2021,07,15, 09, 50, 19), datetime(2021,07,15, 09, 50, 16)]; 
        %datetime(2021,07,15, 09, 52, 28), datetime(2021,07,15, 09, 52, 26);         
        %datetime(2021,07,15, 09, 54, 41), datetime(2021,07,15, 09, 54, 39);
        %];

% Get wav id

nbD =size(dtime,1);
% Path 2 folder
folderIn = ['~/Documents/MPO/BRing/Data/wav/']; 
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/perce/'];
for jj=1:nbD
for ii=1:2
[aziD(jj,ii) detect(jj,ii) ]= fmainBring(arrID{ii},dtime(jj,ii),[folderIn '/' arrID{ii} '/']);
end

locD(jj,:) = beamCrossing(arrID{1},arrID{2}, aziD(jj,1), aziD(jj,2) ,'showFig', 0);
end
%% Show map
pos1 = getArrInfo(arrID{1});
pos2 = getArrInfo(arrID{2});

% Zone de la carte
if strcmp('CLD',arrID{1}) || strcmp('AAV',arrID{1})
lon_min = -65.5;
lon_max = -63.8;
lat_min =  48.9;
lat_max =  49.6;
else    
lon_min = -64.4;
lon_max = -61;
lat_min =  47.5;
lat_max =  48.7;
end

close all
disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);

% figure info
fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen

% Les isobathes à contourer
isobath = [0:10:150];

% Start mapping
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

% Information graphique
c = colorbar('location','eastoutside');
cpos = c.Position;
set(c, 'FontSize', 12, 'FontName','times')
ylabel(c,'Depth(m)')

% Hydrophone location
% ------------ Plot points ----------------
hc1 = m_plot(pos1(2),pos1(1),'ko','MarkerFaceColor','k','MarkerSize',7);
hc2 = m_plot(pos2(2),pos2(1),'ko','MarkerFaceColor','k','MarkerSize',7);

% Plot detection
for ii=1:nbD
    if ~isnan(locD(ii,1))
        m_plot(locD(ii,2),locD(ii,1),'ro','MarkerFaceColor','r','MarkerSize',7)
    end
end

print([getDirectory('fout') '/generalFigures/' 'mappingDectections_' moreInfo  '.png'],'-dpng','-r150')