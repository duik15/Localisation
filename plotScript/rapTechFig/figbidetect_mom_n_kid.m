% This script show the 5 panel diagnostic for a double localisation
% Blue CLD_0950: 19P, 1:56L, 2:28P, 2:44L, 4:41L
% Blue AAV_0950: 16P, 1:56L*, 2:26P, 4:39P
clear all
close all

moreInfo = 'narw_210804T1220_kid';
arrID = {'MLB','PRC'};
%dtime = [ datetime(2021,07,15, 09, 52, 27); datetime(2021,07,15, 09, 52, 25)];
dtime = [ datetime(2021,08,04, 12, 24, 19,30); datetime(2021,08,04, 12, 24, 25)];
%dtime2 = [ datetime(2021,08,04, 12, 24, 18,30); datetime(2021,08,04, 12, 24, 24,30)];

% Parameter
showFig  = [2 3];
imDur = 1.5;
buffer = 0 ;
nbD =1;


% Path 2 folder
%folderWav = {'E:\Bring_dep_2\','F:\Bring_dep_2\'}; 
folderWav = {'/Volumes/BringDD3/Bring_Dep_2/','/Volumes/BringDD4/Bring_Dep_2/'}; 
%folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\perce\']; %spectrogram'];
%folderWav = {'~/Documents/MPO/BRing/Data/wav/', }; 
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/perce/'];

%%
% Parameter
% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [90 210];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 80];           % [dB] C limite pcolor
spgm.im.dur = imDur;%'all';         % [s or 'all'] figure duration
%spgm.im.ovlp = 50;                 % [%] image window overlap
spgm.im.figvision = true ;          % [true false] visiblity of figure before saveas jpf file
% spectrogram window parameters
spgm.win.dur = 2048/10000;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 2;                  % [s] spectrogram zero padding
spgm.im.fmin = spgm.im.freqlims(1);spgm.im.fmax = spgm.im.freqlims(2);

% Parameter for power convertion
convPW.SH =-194;  
convPW.G = 40;
convPW.D = 1;

spgm(2) = spgm(1);
%% Reading file

for ii=1:2
folderIn = [folderWav{ii} arrID{ii} '/'];
% Get file info
[fileListt, wavIDt] = getWavName(dtime(ii), folderIn);
fileList{ii} = fileListt{1}; wavID(ii) = wavIDt;

% Reading wav and extract values
[wav(ii).s,wav(ii).fs, wav(ii).time, wav(ii).info] = readBring([folderIn fileList{ii}], dtime(ii),'duration',spgm(ii).im.dur,'buffer',buffer);
%file_wav = fileList{ifile};     % file name alone
spgm(ii).fs = wav(ii).fs; spgm(ii).ns =  wav(ii).info.Ns; spgm(ii).im.ns= spgm(ii).ns; spgm(ii).im.dur = wav(ii).info.dura;
spgm(ii) = getSpgmWin(spgm(ii));        % Get spectograme windows parameter

wav(ii).pa =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav(ii).s;
end

% Calculate PSD and beamForming
for  ii=1:2

folderIn = [folderWav{ii} arrID{ii} '/'];

% Calculate PSD
[~,freqP{ii},timeP{ii},tmp] = spectrogram(wav(ii).pa(:,1),spgm(ii).win.val,spgm(ii).win.novlp,spgm(ii).win.nfft,spgm(ii).fs);
PSD{ii} = 10*log10(tmp);

% Get angle of arrival
[aziD(ii) detect(ii) ]= fmainBring(arrID{ii},dtime(ii),folderIn);

% Calculate the beamForming
[PSDBF{ii}, timeBF{ii}, freqBF{ii}] = beamForming(arrID{ii}, wav(ii).pa , aziD(ii), spgm(ii),'specmethod','spectro');

end



%% Show spectrogram
if any(showFig == 1)
% Set figure parameter
sp.height=30; sp.width=30; sp.nbx=2; sp.nby=3;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.5; sp.spacey=1;
sp.fracy = [0.2 0.2 0.6];
sp.pos=subplot2(sp);

sp.posMap  = sp.pos{1,3};
sp.posMap(2) = sp.pos{1,3}(2) - 0.05;
sp.posMap(3) = sp.pos{1,3}(3) * 2 + sp.spacex/sp.width;
sp.posMap(4) = sp.pos{1,3}(4) - sp.spacey/sp.height;
      
    
close all

% figure info
kk=1;
fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen


% figure info
kk=1;
%fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
%set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen
fig = figure(1);

for sx=1:2
ax(kk)=axes('position',sp.pos{sx,1},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
pcolor(timeP{sx}, freqP{sx}, PSD{sx}); shading flat;
ylim(spgm(ii).im.freqlims)
set(gca,'FontSize',12)
colormap jet
caxis(spgm(ii).im.clims)

% Label
ax(kk).XTickLabel = [];
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    ylabel('f (Hz)')
end

title({arrID{sx}, 'Hydrophone 1'})

% Panel index
kk=kk+1;
end
%
% Show beamforming
kk=3;
for sx=1:2
ax(kk)=axes('position',sp.pos{sx,2},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
pcolor(timeBF{sx}, freqBF{sx}, squeeze(PSDBF{sx})') ; shading flat;
ylim(spgm(ii).im.freqlims)
set(gca,'FontSize',12)
colormap jet
caxis(spgm(ii).im.clims)

% Label
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    ylabel('f (Hz)')    
end


xlabel('t (s)')
title(['Beam forming at ' num2str(aziD(sx)) '�'])

% Panel index
kk=kk+1;
end

% Colorbar
cbPos = [sp.pos{2,2}(1)+sp.pos{2,2}(3)+0.01 sp.pos{2,2}(2) 0.015 sp.pos{2,2}(4)*2+(sp.spacey/sp.height)];
cb = colorbar;
cb.Position = cbPos;
ax(4).Position = sp.pos{2,2};
ylabel(cb,'PSD (dB re. 1�Pa^2/Hz)','interpreter','tex')    

    
% Show map
kk=5;
ax(kk)=axes('position',sp.posMap,'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');

% And position of location
locD = beamCrossing(arrID{1},arrID{2}, aziD(1), aziD(2) ,'showFig', 0);

% Get position of structure
pos1 = getArrInfo(arrID{1});
pos2 = getArrInfo(arrID{2});

% Zone de la carte
if strcmp('CLD',arrID{1}) || strcmp('AAV',arrID{1})
lon_min = -65.5;
lon_max = -63.5;
lat_min =  49;
lat_max =  50;
rulerXPos=[-65.3 -64.8];
rulerYPos = 49.2;
arrowPos=[-65.3 49.7 1];
else    
lon_min = -64.4;
lon_max = -61;
lat_min =  47.5;
lat_max =  48.7;


end

disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);


% Les isobathes � contourer
isobath = [0:10:150];

% Start mapping
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
%m_ruler;

hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(ax(5),gebco);

% Information graphique
cbPos = [sp.pos{1}(1)+sp.pos{1}(3)+0.01 sp.pos{1}(2) 0.015 sp.pos{1}(4)];
cb2 = colorbar('location','eastoutside');
%cb2Pos = cb2.Position;
%cbPos(4) = cb2Pos(4);
%cb2.Position = cbPos;
%set(gca,'Position',sp.pos{1})
set(cb2, 'FontSize', 12, 'FontName','times')
ylabel(cb2,'Depth(m)')

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

% Ruler
m_ruler([.05 .3],.1,3,'fontsize',8)
%m_northarrow(-65,49.5,1,'type',4,'linewi',.5);

print([getDirectory('fout') '/generalFigures/' 'biDetect5Pan_' moreInfo '.png'],'-dpng','-r150')
end

%% Separeted figures

if any(showFig ==2)
    
% Set figure parameter
clear sp
sp.height=10; sp.width=30; sp.nbx=2; sp.nby=2;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.5; sp.spacey=0.5;
sp.pos=subplot2(sp);
      
    
close all

% figure info
kk=1;
%fig = figure('Units', 'normalized', 'Position', [-1,0,0.8,1]);
%set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen
fig = figure(2);

for sx=1:2
ax(kk)=axes('position',sp.pos{sx,1},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
pcolor(timeP{sx}, freqP{sx}, PSD{sx}); shading flat;
ylim(spgm(ii).im.freqlims)
set(gca,'FontSize',12)
colormap jet
caxis(spgm(ii).im.clims)

% Label
ax(kk).XTickLabel = [];
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    ylabel('$f$ (Hz)','interpreter','latex')
end

title({arrID{sx}, 'Hydrophone 1'})

% Panel index
kk=kk+1;
end
%
% Show beamforming
kk=3;
for sx=1:2
ax(kk)=axes('position',sp.pos{sx,2},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
pcolor(timeBF{sx}, freqBF{sx}, squeeze(PSDBF{sx})') ; shading flat;
ylim(spgm(ii).im.freqlims)
set(gca,'FontSize',12)
colormap jet
caxis(spgm(ii).im.clims)

% Label
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    ylabel('$f$ (Hz)','interpreter','latex')    
end


xlabel('$t$ (s)','interpreter','latex')
title(['Beam forming at ' num2str(aziD(sx)) '�'])

% Panel index
kk=kk+1;
end

% Colorbar
cbPos = [sp.pos{2,2}(1)+sp.pos{2,2}(3)+0.01 sp.pos{2,2}(2) 0.015 sp.pos{2,2}(4)*2+(sp.spacey/sp.height)];
cb = colorbar;
cb.Position = cbPos;
ax(4).Position = sp.pos{2,2};
%ylabel(cb,'PSD (dB re. 1�Pa^2/Hz)','interpreter','tex')    
 ylabel(cb,'PSD (dB re. 1$\mu$ Pa${}^2$/Hz)','interpreter','latex')   
    
%print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetect4Pan_' moreInfo '.png'],'-dpng','-r150')
print([folderOut '/biDetect4Pan_' moreInfo '.png'],'-dpng','-r150')    
end

%% 
if any(showFig==3)
close all
% Set figure parameter
clear sp
clear ax
sp.height=30; sp.width=30; sp.nbx=1; sp.nby=1;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.1; sp.spacey=0.1;
sp.pos=subplot2(sp);

        
% And position of location
[locD distD2]= beamCrossing(arrID{1},arrID{2}, aziD(1), aziD(2) ,'showFig', 0)
distD = mean(distD2);

% Get position of structure
pos1 = getArrInfo(arrID{1});
pos2 = getArrInfo(arrID{2});

% Zone de la carte
if strcmp('CLD',arrID{1}) || strcmp('AAV',arrID{1})
lon_min = -65.5;
lon_max = -64;
lat_min =  49;
lat_max =  50;
else    
    
lon_min = -64.5;
lon_max = -63.5;
lat_min =  48.3;
lat_max =  48.7;
end

disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);


% Les isobathes � contourer
isobath = [0:10:150];

% Start mapping
figure(3)
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
%m_ruler;

hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

% Color bar
cbPos = [sp.pos{1}(1)+sp.pos{1}(3)+0.01 sp.pos{1}(2) 0.015 sp.pos{1}(4)];
cb2 = colorbar('location','eastoutside');
set(cb2, 'FontSize', 12, 'FontName','times')
ylabel(cb2,'Depth(m)')

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

% Ruler
m_ruler([.05 .3],.1,3,'fontsize',8)
%m_northarrow(-65,49.5,1,'type',4,'linewi',.5);

print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetectMap_' moreInfo '.png'],'-dpng','-r150') 
    
end

%% Show approximate pos for kid
if any(showFig ==4)
close all
% Set figure parameter
clear sp
clear ax
sp.height=30; sp.width=30; sp.nbx=1; sp.nby=1;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.1; sp.spacey=0.1;
sp.pos=subplot2(sp);

aziDA = [126 106];
aziDB = [136 113];
        
% And position of location
[locD distD2]= beamCrossing(arrID{1},arrID{2}, aziDA(1), aziDA(2) ,'showFig', 0)
distD = mean(distD2);

% Get position of structure
pos1 = getArrInfo(arrID{1});
pos2 = getArrInfo(arrID{2});

% Limite of map
lon_min = -64.5;
lon_max = -63.8;
lat_min =  48.3;
lat_max =  48.7;

% Find straight line
% Find the pos at X km away
[lon1A,lat1A, phi1A] = m_fdist(pos1(2),pos1(1), aziDA(1), 200 *1000);
[lon2A,lat2A, phi2A] = m_fdist(pos2(2),pos2(1), aziDA(2), 200 *1000);
[lon1B,lat1B, phi1B] = m_fdist(pos1(2),pos1(1), aziDB(1), 200 *1000);
[lon2B,lat2B, phi2B] = m_fdist(pos2(2),pos2(1), aziDB(2), 200 *1000);
[lon2B2,lat2B2, phi2B2] = m_fdist(pos2(2),pos2(1), 109, 200 *1000);
lon1A = lon1A-360; lon1B = lon1B-360; lon2A = lon2A-360; lon2B = lon2B-360;
lon2B2 = lon2B2-360;


lonL1A = linspace(pos1(2),lon1A,1000);
lonL2A = linspace(pos2(2),lon2A,1000);
latL1A = linspace(pos1(1),lat1A,1000);
latL2A = linspace(pos2(1),lat2A,1000);
lonL1B = linspace(pos1(2),lon1B,1000);
lonL2B = linspace(pos2(2),lon2B,1000);
latL1B = linspace(pos1(1),lat1B,1000);
latL2B = linspace(pos2(1),lat2B,1000);
latL2B2 = linspace(pos2(1),lat2B2,1000);
lonL2B2 = linspace(pos2(1),lon2B2,1000);

[lonA,latA] = polyxpoly2(lonL1A,latL1A,lonL2A,latL2A);
[lonB,latB] = polyxpoly2(lonL1B,latL1B,lonL2B,latL2B);
[lonB2,latB2] = polyxpoly2(lonL1B,latL1B,lonL2B2,latL2B2);
[lonm,latm] = polyxpoly2(lonL1B,latL1B,lonL2A,latL2A);


% Calcul the B position
dtW = 2.3;
dW = dtW * 1.470;

for ii=1:length(lonL1B)
distA = m_idist(pos1(2),pos1(1),lonA,latA) / 1000;
dist1 = m_idist(pos1(2),pos1(1),lonL1B(ii),latL1B(ii)) / 1000  ;
dist2 = m_idist(lonA,latA,lonL1B(ii),latL1B(ii)) / 1000;
sumD = (dist1 + dist2 );
diffD(ii) = sumD - distA;
end

[valm, indm] = min(abs(diffD) ); 
[val, indB] = min(abs(diffD -  dW));
[vals inds] = sort(abs(diffD -  dW));
%
disp('Start mapping')
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);


% Les isobathes � contourer
isobath = [0:10:150];

% Start mapping
figure(3)
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
%m_ruler;

hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

% Color bar
cbPos = [sp.pos{1}(1)+sp.pos{1}(3)+0.01 sp.pos{1}(2) 0.015 sp.pos{1}(4)];
cb2 = colorbar('location','eastoutside');
set(cb2, 'FontSize', 12, 'FontName','times')
ylabel(cb2,'Depth(m)')

% Hydrophone location
% ------------ Plot points ----------------
hc1 = m_plot(pos1(2),pos1(1),'ko','MarkerFaceColor','k','MarkerSize',7);
hc2 = m_plot(pos2(2),pos2(1),'ko','MarkerFaceColor','k','MarkerSize',7);

% Plot detection
%for ii=1:nbD
%    if ~isnan(locD(ii,1))
%        m_plot(locD(ii,2),locD(ii,1),'ro','MarkerFaceColor','r','MarkerSize',7)
%    end
%end
m_plot(lonA,latA,'ro','MarkerFaceColor','r','MarkerSize',7)
%m_plot(lonB,latB,'ro','MarkerFaceColor','r','MarkerSize',7)
!m_plot(lonL1B(inds(1)),latL1B(inds(1)),'co','MarkerFaceColor','r','MarkerSize',7)
m_plot(lonL1B(inds(2)),latL1B(inds(2)),'go','MarkerFaceColor','g','MarkerSize',7)
m_plot(lon2B2,lat2B2,'go','MarkerFaceColor','g','MarkerSize',7)
%m_plot(lonL1B(inds(2)),latL1B(inds(2)),'o','color',[66 167 1]/255,'MarkerFaceColor',[66 167 1]/255,'MarkerSize',7)


% Line of detection
ii=1;
hl1= m_plot(lonL1A(ii,:), latL1A(ii,:),'-.','color',0.3*[1 1 1],'lineWidth',1.5);
hl1= m_plot(lonL2A(ii,:), latL2A(ii,:),'-.','color',0.3*[1 1 1],'lineWidth',1.5);
hl1= m_plot(lonL1B(ii,:), latL1B(ii,:),'--','color',0.3*[1 1 1],'lineWidth',1.5);
%hl1= m_plot([lonL1B(inds(2)) lonm ],[latL1B(inds(2)) latm],'-.','color','g','lineWidth',1.5);
%hl1= m_plot(lonL2B(ii,:), latL2B(ii,:),'-.','color',0.3*[1 1 1],'lineWidth',1.5);



% Ruler
m_ruler([.05 .3],.1,3,'fontsize',8)
%m_northarrow(-65,49.5,1,'type',4,'linewi',.5);

print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetectMap_analyzed_' moreInfo '.png'],'-dpng','-r150') 
    
%%
% Segment 
hl1= m_plot([lonL1B(inds(2)) pos1(2)],[latL1B(inds(2)) pos1(1)],'-.','color',[0.4940 0.1840 0.5560],'lineWidth',1.5);
hl1= m_plot([lonL1B(inds(2)) lonA],[latL1B(inds(2)) latA],'-.','color',[0.4940 0.1840 0.5560],'lineWidth',1.5);
hl1= m_plot([pos1(2) lonA],[pos1(1) latA],'-.','color',[0.8500 0.3250 0.0980],'lineWidth',1.5);
%hl1= m_plot([lonL1B(inds(2)) lonm ],[latL1B(inds(2)) latm],'-','color','g','lineWidth',2.5);
print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetectMap_analyzedLine_' moreInfo '.png'],'-dpng','-r150') 

end