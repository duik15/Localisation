%function [loc bc] = beamCrossing(arrID1,arrID2,a1,a2,varargin)
% This script make a triangularization of two beam to find back to position
%[loc] = beamCrossing(pos1,pos2,a1,a2,varargin)
%https://www.mathworks.com/help/map/ref/polyxpoly.html
% Need update
% need to add incertainty
% need to be able to send 2 datetime, find best angle for each and than beamcrossing
% Last moif @ kevDuquette  3 nov 2021


%debug:
arrID1 = 'MLB';
arrID2 = 'PRC';
a1=100;
a2=70;

showFig = [1];
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/honguedo/'];%['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\honguedo\'];
outName = [arrID1 '_' num2str(a1) 'deg_' arrID2 '_' num2str(a2) 'deg' ];
mkdir2(folderOut)

pos1 = getArrInfo(arrID1);
pos2 = getArrInfo(arrID2);

% Default parameter
limType = 'fixe';
inc=2;

%% Varagin
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'cellnumber'
            cellNumber = varargin{2};
        case 'showfig'
            showFig = varargin{2};
        case 'limtype'
            limType = varargin{2};
            
        otherwise
            error(['Can''t understand property: ' varargin{1}])
    end
    varargin(1:2)=[];
end


%% Find angle curve

[lon1m,lat1m, phi1] = m_fdist(pos1(2),pos1(1), a1, 100 *1000);
[lon2m,lat2m, phi2] = m_fdist(pos2(2),pos2(1), a2, 100 *1000);

lon1m = lon1m - 360;
lon2m = lon2m -360;


lon1 = linspace(pos1(2),lon1m,1000);
lon2 = linspace(pos2(2),lon2m,1000);
lat1 = linspace(pos1(1),lat1m,1000);
lat2 = linspace(pos2(1),lat2m,1000);

[loni,lati] = polyxpoly2(lon1,lat1,lon2,lat2);

loc = [lati, loni];

%% Calcul incertitude zone
da = [-inc inc];
ii=0;
for i1 = 1:2
    for i2 = 1:2
        ii= ii+1;
        [lon1m,lat1m, phi1] = m_fdist(pos1(2),pos1(1), a1 + da(i1), 100 *1000);
        [lon2m,lat2m, phi2] = m_fdist(pos2(2),pos2(1), a2 + da(i2), 100 *1000);
        
        lon1m = lon1m - 360;
        lon2m = lon2m -360;
        
        lon11 = linspace(pos1(2),lon1m,1000);
        lon22 = linspace(pos2(2),lon2m,1000);
        lat11 = linspace(pos1(1),lat1m,1000);
        lat22 = linspace(pos2(1),lat2m,1000);

        [londa(ii),latda(ii)] = polyxpoly2(lon11,lat11,lon22,lat22);        
        
        % SAving value for plot
        if i2==1
            lon1da(i1,:) = lon11; lat1da(i1,:) = lat11;
        end
        if i1==1
            lon2da(i2,:) = lon22; lat2da(i2,:) = lat22;
        end
        
    end
end

% Change order so it turn in circle
londa = londa([1 2 4 3]);
latda = latda([1 2 4 3]);

%% Find the limit of the zone to plot

% Zone de la carte / Limits
if strcmp(limType , 'fixe')
    winSize = 0.3;
    if strcmp(arrID1,'PRC') || strcmp(arrID2,'PRC')
        % Zone de la carte
        lon_min = -64.4;
        lon_max = -63.9;
        lat_min =  48.4;
        lat_max =  48.7;
    else strcmp(arrID1,'AAV') || strcmp(arrID2,'AAV')
        
    end
    
else strcmp(limType,'auto')
    %lon_min = -max(abs([pos1(:,2) pos2(:,2)])) - winSize;
    %lon_max = -min(abs([pos1(:,2) pos2(:,2)])) + winSize ;
    %lat_min =  min([pos1(:,1) pos2(:,1)]) -  winSize/2;
    %lat_max =  max([pos1(:,1) pos2(:,1)]) + winSize/2;
end

%{
    dlon = abs(diff([pos1(2) pos2(2)]));
    dlat = abs(diff([pos1(1) pos2(1)]));
    dl = max([dlon dlat])
    lon_min = min([pos1(2) pos2(2)]) - dlon/2;
    lon_max = max([pos1(2) pos2(2)]) + dlon/2;
    lat_min = min([pos1(1) pos2(1)]) - dlat/4;
    lat_max = max([pos1(1) pos2(1)]) + dlat/4;
%}


%% Figure
%if any(showFig ==1)
    close all
    disp('Start mapping')
    % Load map
    gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
    gebco = load([getDirectory('map') 'gebco_colormap.dat']);

    
    % Les isobathes � contourer
    isobath = [0:5:150];
    
    fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
    set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen
    
    m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
    m_grid('fontsize',12,'linestyle','none');
    
   
    
    hold on
    m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
    
    
    % Information graphique
    c = colorbar('location','eastoutside');
    cpos = c.Position;
    set(c, 'FontSize', 12, 'FontName','times')
    ylabel(c,'Depth(m)')
    
    % Un beau colormap pour la bathymetrie qui vient de GEBCO
    colormap(gebco);
    
    %xlabel('Longitude', 'FontSize', 12, 'FontName','times')
    %ylabel('Latitude', 'FontSize', 12, 'FontName','times')
    
    m_gshhs_f('patch',[.7 .7 .7]); % coastlines
    

    % Hydrophone location
    % ------------ Plot points ----------------
    hc1 = m_plot(pos1(2),pos1(1),'ro','MarkerFaceColor','r','MarkerSize',7);
    hc2 = m_plot(pos2(2),pos2(1),'ro','MarkerFaceColor','r','MarkerSize',7);
    
    % Plot line crossing
    %l1 = m_plot(lon1,lat1,'-.k');
    %l2 = m_plot(lon2,lat2,'-.k');
    
 
   % Cone
   for ii=1:2
    hl1= m_plot(lon1da(ii,:), lat1da(ii,:),'-.','color',0.3*[1 1 1],'lineWidth',1.5)
    hl2= m_plot(lon2da(ii,:), lat2da(ii,:),'-.','color',0.3*[1 1 1],'lineWidth',1.5)
   end
    % Plot incertitude
    %for ii=1:4
    %    pda = m_plot(londa,latda,'ko','MarkerFaceColor','k','MarkerSize',4);
    %end
    % Plot zone
    m_line([londa londa(1)],[latda latda(1)],'color',[0 0 0])
   % m_patch([londa(1,:) londa(1,:)],[latda(1,:) latda(1,:)],0.3)
      
    % Plot the retrieve point
    pp = m_plot(loni,lati,'ko','MarkerFaceColor','k','MarkerSize',5);
    
    %%
    % Legend
    %leg = m_legend([hc1,pp], 'Array1','Retreive signal')
    %leg.Position = [49,-64,10,4]
    %leg.Position = [49,5,25,3]
    
    set(gca,'YMinorTick','on','XMinorTick','on')
    print('-dpng', '-r300',[folderOut 'beamCrossing_' outName '.png']);
%end
