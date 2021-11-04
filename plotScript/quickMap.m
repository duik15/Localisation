function quickMap(lat,lon, varargin)
% This function quicly plot positions
% quickMap(lat,lon,'type')
% typre of map: online / geoshow / mmap
disp('Printing map....')


type = 'mmap';
typePlot = 'point' ;  

%% Varagin
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                time = varargin{2};
            case 'type'
                type = varargin{2};
            case 'typeplot'
                 typePlot = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end




% --------------------- Quick GeoShow ----------------------
if strcmp(type , 'geoshow')
    
% Plot limiit
[latlim, lonlim] = geoquadline(lat, lon);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, .01, .01);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, .1, .1);

% Load coasyline
load coastlines

fig = figure;
pos = fig.Position;
fig.Position = [300 300 1.25*pos(3) 1.25*pos(4)];
ax = usamap(latlim, lonlim);
setm(ax, 'MLabelParallel', 43.5)

geoshow(coastlat,coastlon, ...
    "DisplayType","polygon","FaceColor",[0.45 0.60 0.30])

geoshow(lat, lon, 'DisplayType' , 'point')
end


%  ------------------ Online ----------------------
if strcmp(type , 'online')
    webmap('openstreetmap')
    wmline(lat,lon, 'Color', 'red')
end


% ------------- MMAP -----------------------------
if strcmp(type , 'mmap')
    
    try
        % Read bathymetry
        load Documents/MATLAB/Oceanographie/m_map/mmap_ex/GSL_bathy_500m.mat;
        load Documents/MATLAB/Oceanographie/m_map/mmap_ex/gebco_colormap.dat;
    catch
        try
            % Read bathymetry
            load GSL_bathy_500m.mat;
            load gebco_colormap.dat;
        catch
            error("Cant not load bathy file.")
        end
    end
    
    
% Default parameter
winSize = 0.05;

% Zone de la carte
lon_min = min(lon)-winSize;
lon_max = max(lon)+winSize;
lat_min =  min(lat)-winSize;
lat_max =  max(lat)+winSize;

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

% ------------ Plot points ----------------

for i=1:length(lon)
hr=   m_plot(lon(i),lat(i),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','k','MarkerSize',4) ;
end
end

%print('-dpng', '-r300','figures/WBR_map_honguedo_transNCirc.png');    
disp('Done!')
end
