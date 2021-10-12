% This script takes the ping input and boat input to make a angle of
% arrival line
clear all

arrID = 'AAV';

% Loading data
folderPing = '/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/OffSet23_voc3_Ns2e14_150-200hz_AvCircle/';
p = load([folderPing 'dataAngleOfArrival_OffSet23_voc3_Ns2e14_150-200hz_AvCircle.mat'])

b2 = load('/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/OffSet23_voc3_Ns2e14_150-200hz_boat1minute/dataAngleOfArrival_OffSet23_voc3_Ns2e14_150-200hz_boat1minute.mat');
folderBoat = '/Users/Administrator/Documents/MPO/BRing/Data/gpxTrack/';
b = load([folderBoat 'aavCircleTrack.mat']);
b.angleM = b2.angleM;
b.angleA = b2.angleA;

% Get the array location
arrLoc = getArrLoc(arrID);

% Find the ping angle
for  i =1:length(p.ptime)
[distp(i),apcp(i),appc(i)] = m_idist(arrLoc(2),arrLoc(1),-p.ploc(i,2),p.ploc(i,1));
end
p.angleLL = 360 -apcp;

% Find the angle of the boat sound
for  i =1:length(b.time)
[distb(i),abcp(i),abpc(i)] = m_idist(arrLoc(2),arrLoc(1),b.lon(i),b.lat(i));
end
b.angleLL = 360 - abcp; 



%% Map figure

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
lon_min = min(b.lon)-winSize;
lon_max = max(b.lon)+winSize;
lat_min =  min(b.lat)-winSize;
lat_max =  max(b.lat)+winSize;

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

ha=   m_plot(arrLoc(2),arrLoc(1),'*','MarkerEdgeColor','k',...
    'MarkerFaceColor','k','MarkerSize',6) ;

hold on

% Boat location
hb=   m_plot(b.lon(:),b.lat(:),'o','MarkerEdgeColor','b',...
    'MarkerFaceColor','b','MarkerSize',4) ;

%for i=1:length(b.lon)
%m_text(b.lon(i)+ 0.0001 ,b.lat(i)+ 0.0001,num2str(i),'color','b');
%end

% Ping location
hp=   m_plot(-p.ploc(:,2),p.ploc(:,1),'+','MarkerEdgeColor','r',...
    'MarkerFaceColor','r','MarkerSize',4) ;

%for i=1:length(p.ptime)
%m_text(-p.ploc(i,2)+ 0.0001 ,p.ploc(i,1)+ 0.0001,num2str(i),'color','r');
%end

legend([ha,hb,hp],'Center Array','Ping Location','Boat Location')
print('-dpng', '-r150',[folderPing 'pingNBoatMap.png']);    
disp('Done!')



%% Figure

nbP = length(p.ptime);
nbB = length(b.time);

if sum(diff(p.angleM(3:end-2)))< 0
    pNb = nbP:-1:1;
    angleA = fliplr(aMax3);
else
    pNb = 1:nbP;
end    
    
% Linéaire regression
xp = 1:360;
ppfit = polyfit(360 - p.angleLL,p.angleM,1);
angleFitp = polyval(ppfit,xp);  

xb = 1:360;
pbfit = polyfit(360 - b.angleLL,b.angleM,1);
angleFitb = polyval(pbfit,xb);  

figure(1)
ax = gca;
colorO = ax.ColorOrder;
hp = plot(360 -p.angleLL,p.angleM,'-.o','color',colorO(1,:));
hold on
hb = plot(360 - b.angleLL,b.angleM,'o','color',colorO(2,:));
hb2 = plot(360 - b.angleLL,b.angleA(:,2),'*','color',colorO(2,:));
hb3 = plot(360 - b.angleLL,b.angleA(:,3),'+','color',colorO(2,:));
%hpTh = plot(x,0:36:360-36,'k-');
hbFit = plot(xb,angleFitp,'color',colorO(1,:));
hpFit = plot(xb,angleFitb,'color',colorO(2,:));

% Plot some side lobe
%for ii=2:2
%   hSl = plot(x,angleA(:,ii),'o') 
%end

xlabel('Theoritical angle')
ylabel('Beamforming angle')
%ax.XTickLabel = pNb;
leg= legend([hp,hb,hb2, hb3, hpFit,hbFit],{'Ping 1st lobe','Boat 1st lobe','Boat 2nd lobe','Boat 3rd lobe','Ping lin. fit','Boat lin. fit'},'Location','northwest')

print([folderPing 'angleOfArrival_wBoat.png'],'-r150','-dpng')
