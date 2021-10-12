% Showing the hydrophone location on a map

arrID = 'AAV';

switch lower(arrID)
    case 'aav'
        [arrLoc offset]= getArrLoc(arrID);
end
        
lim = 12;
Nc = 20;
R = 10;
theta = 0: 2*pi/(Nc) : 2*pi-2*pi/(Nc);

% Version x,y
xc = R*sin(theta +  offset * pi /180);
yc = R*cos(theta +  offset * pi /180);


%% Convert to lat, lon
for i=1:Nc
[hlon2(i,1), hlat(i,1), a21(i,1)] = m_fdist(-arrLoc(2), arrLoc(1),theta(i) * 180  / pi+  offset -90,R);        
[hlont2(i,1), hlatt(i,1), at21(i,1)] = m_fdist(-arrLoc(2), arrLoc(1),theta(i) * 180  / pi+  offset -90,R+1);        
end
hlon = -hlon2;
hlont = -hlont2;

%% Plot in x,y
figure(1)

h1 = plot(xc,yc,'-ok');
hold on
h2 = plot([xc(1) xc(end)],[yc(1) yc(end)],'-ok');

ax=gca;

dx = 0.2;
for ic = 1:length(xc)
text(xc(ic)+dx,yc(ic)+dx,num2str(ic))
end

% Plot arrow
hq = quiver(0,-lim,0,2*lim+2.5,'-k');
hq = quiver(-lim,0,2*lim+2.5,0,'-k');

% Plot N,E,S,W
text(0.7,11,'N','fontsize',14)
text(0.7,-11,'S','fontsize',14)
text(-11,0.7,'W','fontsize',14)
text(11,0.7,'E','fontsize',14)

% Figure parameter
xlim([-lim, lim])
ylim([-lim, lim])
axis square
xlabel('X (m)')
ylabel('Y (m)')
grid on
box on
print('-dpng','-r150',['~/Documents/MPO/BRing/Data/results/generalFigures/' arrID '_hydroMap.png '])

%% Plot lat,lon 1

%xliml = arrLoc(1) * [1 1] + 0.0001 * [-1 1];
%yliml = arrLoc(2) * [1 1] +  0.0001 * [-1 1] ;

% Plot limiit
[latlim, lonlim] = geoquadline(hlat, hlon);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, 0.00005 , .00005);

clf
figure(1)

h1 = plot(hlat,hlon,'-ok');
hold on
h2 = plot([hlat(1) hlat(end)],[hlon(1) hlon(end)],'-ok');

ax=gca;

% Plot the hydrophone number

for i=1:Nc
text(hlat(i),hlont(i), num2str(i))
end
%
%
% Plot arrow
hq = quiver( arrLoc(1), lonlim(1), 0, diff(lonlim)+0.00003,'-k');
hq = quiver( latlim(1), arrLoc(2), diff(latlim)+0.00003, 0,'-k');
%hq = quiver( ax.XLim(1), sum(ax.YLim)/2,  diff(ax.XLim)+0.00005 ,0,'-k');
%hq = quiver(-lim,0,2*lim+2.5,0,'-k');
%
% Plot N,E,S,W
%text(arrLoc(1) , latlim(2),'N','fontsize',14)
%text(sum(ax.XLim)/2 ,ax.YLim(1)-0.00004,'S','fontsize',14)


% Figure parameter
xlim( latlim )
ylim( lonlim )
axis square
xlabel('Latitude ')
ylabel('Longitude')
grid on
box on
print('-dpng','-r150',['~/Documents/MPO/BRing/Data/results/generalFigures/' arrID '_hydroMap_latlon1.png '])


%% Plot lat,lon 1

%xliml = arrLoc(1) * [1 1] + 0.0001 * [-1 1];
%yliml = arrLoc(2) * [1 1] +  0.0001 * [-1 1] ;

% Plot limiit
[latlim, lonlim] = geoquadline(hlat, hlon);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, 0.00005 , .00005);

clf
figure(1)
for i=1:Nc
text(hlat(i),hlont(i), num2str(i))
end

ax=gca;


% Figure parameter
xlim( latlim )
ylim( lonlim )
axis square
xlabel('Latitude ')
ylabel('Longitude ')
grid on
box on

print('-dpng','-r150',['~/Documents/MPO/BRing/Data/results/generalFigures/' arrID '_hydroMap_latlon_number.png '])




%% Plot in lat,lon

% Plot limiit
[latlim, lonlim] = geoquadline(hlat, hlon);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, 0.00005 , .00005);
%[latlim, lonlim] = geoquadline(route.Latitude, route.Longitude);
%[latlim, lonlim] = bufgeoquad(latlim, lonlim, .05, .05);


%
fig = figure;
pos = fig.Position;
%fig.Position = [300 300 1.25*pos(3) 1.25*pos(4)];
ax = usamap(latlim, lonlim);
setm(ax, 'MLabelParallel', 43.5)
geoshow(hlat', hlon', 'DisplayType' , 'point')

% Tick
ax.XTick = [hlat(10) hlat(16)]

xlabel('Latitude')
xlabel('Longitude')
northarrow('latitude', latlim(2)-0.00005,'longitude',lonlim(2)-0.00005)

print('-dpng','-r150',['~/Documents/MPO/BRing/Data/results/generalFigures/' arrID '_hydroMap_latlon.png '])

 
%% Mmap

    
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
lon_min = min(hlon)-winSize;
lon_max = max(hlon)+winSize;
lat_min =  min(hlat)-winSize;
lat_max =  max(hlat)+winSize;

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

% ------------ Plot ctd ----------------
% Perpendicular transect
for i=1:length(lon)
hr=   m_plot(hlon(i),hlat(i),'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','k','MarkerSize',4) ;
end


%print('-dpng', '-r300','figures/WBR_map_honguedo_transNCirc.png');    
disp('Done!')

