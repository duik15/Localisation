function [time, loc, lat, lon, angle, dist] = loadGPXTrack(timeIn, gID,varargin)
% [gtime gloc] = loadGPXTrack(gID,varargin)
% [gtime gloc] = loadGPXTrack(gID,'istart',istart,'iend',iend,'type','online')
% Varargin option:
%   type : online = open street map show
%          geoshow= quick matlab blanck map
%          mmap   = heavy and full resolution mmap
%   gID :   aav, cld, mlb, prc -> Track around the array
%            file -> loadGPXTrack('file',file2load)

% Debug
%gID  = 'aav';
%timeIn = [];

% defautl parameter
%nbPoint = 10;
showFig = 0;
printFig = [];
saveTrack = true;

% Get the path to data
compName = getMachine();

switch lower(compName)
    case 'mac'
        path2track = '/Users/Administrator/Documents/MPO/BRing/Data/boatTrack/';
    case 'dellkev'
        path2track = 'C:\Users\duquettek\Documents\BRing\Data\boatTrack\';
    otherwise
        warning('Path to boat track is not fund. Please make sure to specify the path in the varagin option or add your computer to the list.')
end

%%
switch lower(gID)
    case 'aav'
        fileName = 'aavCircle.gpx';
        istart    = 1600;
        iend      = 1900;
    case 'cld'
        fileName = 'cldCircleTrack.gpx';
        istart    = 1;
        iend      = 300;
    case 'mlb'
        fileName = 'mlbPrcTrack.gpx';
        istart    = 110;
        iend      = 1350;
    case 'prc'
        fileName= 'mlbPrcTrack.gpx';
        istart    = 650;
        iend      = 850;
        
    case 'file'
        file2load = varargin{1};
        varargin(1)=[];
end

while ~isempty(varargin)
        switch lower(varargin{1})
            case 'nbP'
                nbP = varargin{2};
            case 'istart'
                istart = varargin{2};
            case 'iend'
                iend = varargin{2};
            case 'showfig'
                showFig = varargin{2};
            case 'printfog'
                printFig = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end

%%

if ~exist('file2load')
    file2load = [path2track fileName]
end

% Load the route     
routeF = gpxread(file2load);

if ~exist('istart')
    istart  = 1;
end
if ~exist('iend')
    istart  = length(routeF);
end

% Selected index
route.lat = routeF.Latitude(istart:iend);
route.lon = routeF.Longitude(istart:iend);

% Date and Time
timeStr = strrep(routeF.Time(istart:iend),'Z','');
dateN = datenum(timeStr,'yyyy-mm-ddThh:MM:ss');
route.Time = datetime(dateN,'ConvertFrom', 'datenum');

% Find the time asked
if isempty(timeIn)
    index= 1:length(route.Time);
else
    index = nan(numel(timeIn),1);
    
    for i=1:numel(timeIn)
        [~,index(i)] = min(abs(route.Time - timeIn(i))); % Index of the depth
    end
end
%index = unique(index);

% Saving value
time = route.Time(index);
lat = route.lat(index);
lon = route.lon(index);
loc = [lat',  lon'];

[angle dist]= getRealAngle(gID, lat,lon);

if saveTrack == true
%save(['/Users/Administrator/Documents/MPO/BRing/Data/gpxTrack/aavCircleTrack.mat'],'route')
save(['/Users/Administrator/Documents/MPO/BRing/Data/boatTrack/' gID 'CircleTrack.mat'],'lat','lon','time','angle','dist')
end 

%----------------------  Figure  ------------------------------
if any(showFig ==1 )
    
% Plot limiit
[latlim, lonlim] = geoquadline(route.lat, route.lon);
[latlim, lonlim] = bufgeoquad(latlim, lonlim, .01, .01);    
    
fig = figure;
pos = fig.Position;
fig.Position = [300 300 1.25*pos(3) 1.25*pos(4)];
ax = usamap(latlim, lonlim);
setm(ax, 'MLabelParallel', 43.5)
geoshow(glat, glon)


if ~isempty(printFig)
print('-dpng','-r150',[printFig '/trackMap.png'])    
%wmmarker(wpt, 'FeatureName', wpt.Name, 'OverlayName', 'Waypoints')
end 
end 


if any(showFig  ==2)

    
disp('Printing map')
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
winSize = 0.01;

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
  

% ------------ Plot position ----------------

%ha=   m_plot(arrLoc(2),arrLoc(1),'*','MarkerEdgeColor','k',...
%    'MarkerFaceColor','k','MarkerSize',6) ;

hold on

% Boat location
hb=   m_plot(lon,lat,'o','MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',4) ;

for i=1:length(lon)
m_text(lon(i)+ 0.0001 ,lat(i)+ 0.0001,num2str(i),'color','b');
end

legend(hb,'Boat Location')
%print('-dpng', '-r300','figures/WBR_map_honguedo_transNCirc.png');    
disp('Done!')



    
end

%end

