function s = loadGPXTrack(timeIn, gID,varargin)
%
% s = loadGPXTrack(gID,varargin)
% s = loadGPXTrack(gID,'istart',istart,'iend',iend,'type','online')
% Create a new .mat 
% ~ = loadGPXTrack([],'AAV','saveTrack',true,'istart',1,'iend',0)
% Varargin option:
%   type : online = open street map show
%          geoshow= quick matlab blanck map
%          mmap   = heavy and full resolution mmap
%   gID :   aav, cld, mlb, prc -> Track around the array
%            file -> loadGPXTrack('file',file2load)
% Old :[time, loc, lat, lon, angle, dist] = loadGPXTrack(timeIn, gID,varargin)
% Debug
%gID  = 'aav';
%timeIn = [];

% defautl parameter
%nbPoint = 10;
showFig = [];
printFig = false;
saveTrack = false;
gpxToolBox = false;
moreInfo = '';

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

switch lower(gID)
    case 'aav'
        if gpxToolBox == true
            fileName = 'aavCircleTrack.gpx';
        else
            fileName = 'aavFullTrack.mat';
        end
        istart    = 1600;
        iend      = 1900;
    case 'cld'
        if gpxToolBox == true
            fileName = 'cldCircleTrack.gpx';
        else
            fileName = 'cldFullTrack.mat';
        end
        istart    = 1;
        iend      = 300;
    case 'mlb'
        if gpxToolBox == true
            fileName = 'mlbPrcTrack.gpx';
        else
            fileName = 'mlbFullTrack.mat';
        end
        
        istart    = 110;
        iend      = 350;
    case 'prc'
        if gpxToolBox == true
            fileName = 'mlbPrcTrack.gpx';
        else
            fileName = 'prcFullTrack.mat';
        end
        
        istart    = 650;
        iend      = 850;
        
    case 'file'
        file2load = varargin{1};
        varargin(1)=[];
end
%%
while ~isempty(varargin)
        switch lower(varargin{1})
            case 'nbP'
                nbP = varargin{2};
            case 'istart'
                istart = varargin{2};
            case 'iend'
                iend = varargin{2};
            case 'timemin'
                timemin = varargin{2};
                istart = 1;
            case 'timemax'
                timemax = varargin{2};
                iend=0;
            case 'showfig'
                showFig = varargin{2};
            case 'printfig'
                printFig = varargin{2};
            case 'savetrack'
                saveTrack = varargin{2};
            case 'path2track'
                path2track = varargin{2};
            case 'folderout'
                outPath = varargin{2};
            case 'gpxToolBox'
                gpxToolBox = varargin{2};
            case 'moreInfo'
                moreInfo = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end

%% 

if ~exist('file2load')
    file2load = [path2track fileName];
end

% Load the route 
if gpxToolBox == true
    routeF = gpxread(file2load);
else
    routeF = load(file2load);
end

if ~exist('istart')
    istart  = 1;
end
if ~exist('iend') || iend ==0
    if gpxToolBox == true
        iend  = length(routeF.Latitude);
    else
        iend = length(routeF.lat);
    end
    
end

if gpxToolBox == true
% Selected index
route.lat = routeF.Latitude(istart:iend);
route.lon = routeF.Longitude(istart:iend);

% Date and Time
timeStr = strrep(routeF.Time(istart:iend),'Z','');
dateN = datenum(timeStr,'yyyy-mm-ddThh:MM:ss');
route.Time = datetime(dateN,'ConvertFrom', 'datenum');
else
   route.Time  = routeF.time(istart:iend);
   route.lat  = routeF.lat(istart:iend);
   route.lon  = routeF.lon(istart:iend);
    
end
% Find the time asked
if isempty(timeIn)
    index= 1:length(route.Time);
else
    index = nan(numel(timeIn),1);
    
    for i=1:numel(timeIn)
        [~,index(i)] = min(abs(route.Time - timeIn(i))); % Index of the depth
    end
end
index = unique(index,'stable');

% Add another time filter if timemax and timemin is specify
if exist('timemin')
    indMin = min(find(route.Time > timemin));
    index = index(index> indMin);
end
if exist('timemax')
    indMax = max(find(route.Time < timemax));
    index = index(index < indMax);
end

% Saving value
time = route.Time(index);
lat = route.lat(index)';
lon = route.lon(index)';
loc = [lat,  lon];

% Get the real position
[angle, dist]= getRealAngle(gID, lat,lon);

% Create a structure
s.time = time; s.lat = lat; s.lon = lon; s.loc=loc; s.angle = angle; s.dist = dist;

if saveTrack == true
%save(['/Users/Administrator/Documents/MPO/BRing/Data/gpxTrack/aavCircleTrack.mat'],'route')
save([path2track gID 'Track.mat'],'lat','lon','time','angle','dist')
%save([path2track gID 'FullTrack.mat'],'lat','lon','time','angle','dist')
disp('Done saving track!')
end 

%% ----------------------  Figure  ------------------------------

if printFig == true
   if ~exist('outPath')
       outPath = [path2track 'trackMap\']; 
   end
   if ~isfolder(outPath); disp(['Creating output folder: ' outPath]); mkdir(outPath); end
end


if any(showFig ==1)
    % Plot limit
    %dlat = 0.00005;%max(lat) - min(lat);
    %dlon = 0.00005;%max(abs(lon)) - min(abs(lon));
    figure(1)    
    h1 = plot(lon,lat,'-ok');
    ax=gca;
    grid on;
    xlabel('Longitude')
    ylabel('Latitude')
    if printFig == true
        print('-dpng','-r150',[outPath gID 'Track.png'])  
    end
end

if any(showFig==2)

    figure(1)   
    subplot(3,2,1:2)
    [ax1,h1,h2] = plotyy(time, dist,time,angle);
    xlabel('Time')
    ylabel(ax1(1),'Distance(m)')
    ylabel(ax1(2),'Azimute($^\circ$)','interpreter','latex')
    
    subplot(3,2,3:6)
    h1 = plot(lon,lat,'-ok');
    ax=gca;
    grid on;
    xlabel('Longitude')
    ylabel('Latitude')
    axis square
    if printFig == true
        print('-dpng','-r150',[outPath gID 'TrackNdistance.png'])  
    end
    
end



if any(showFig ==3 )
        
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


if any(showFig  ==4)

    
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

