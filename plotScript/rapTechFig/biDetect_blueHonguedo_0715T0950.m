% This script show the 5 panel diagnostic for a double localisation
% Blue CLD_0950: 19P, 1:56L, 2:28P, 2:44L, 4:41L
% Blue AAV_0950: 16P, 1:56L*, 2:26P, 4:39P
clear all
close all

moreInfo = 'blue150721';
arrID = {'CLD','AAV'};  
dtime = [ datetime(2021,07,15, 09, 52, 27); datetime(2021,07,15, 09, 52, 25)];
%dtime = [ datetime(2021,07,15, 09, 50, 17); datetime(2021,07,15, 09, 50, 15)];
%dtime = [ datetime(2021,07,15, 09, 52, 17); datetime(2021,07,15, 09, 50, 15)];

% Parameter
showFig  = [1 2 3];
imDur = 3;
buffer = 0 ;
nbD =1;


% Path 2 folder
folderWav = ['~/Documents/MPO/BRing/Data/wav/']; 
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/honguedo/'];

% Parameter
% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [0 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 100];           % [dB] C limite pcolor
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
% Reading file

for ii=1:2
folderIn = [folderWav arrID{ii} '/'];
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

folderIn = [folderWav arrID{ii} '/'];

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
title(['Beam forming at ' num2str(aziD(sx)) '°'])

% Panel index
kk=kk+1;
end

% Colorbar
cbPos = [sp.pos{2,2}(1)+sp.pos{2,2}(3)+0.01 sp.pos{2,2}(2) 0.015 sp.pos{2,2}(4)*2+(sp.spacey/sp.height)];
cb = colorbar;
cb.Position = cbPos;
ax(4).Position = sp.pos{2,2};
ylabel(cb,'PSD (dB re. 1µPa^2/Hz)','interpreter','tex')    

    
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


% Les isobathes à contourer
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
title(['Beam forming at ' num2str(aziD(sx)) '°'])

% Panel index
kk=kk+1;
end

% Colorbar
cbPos = [sp.pos{2,2}(1)+sp.pos{2,2}(3)+0.01 sp.pos{2,2}(2) 0.015 sp.pos{2,2}(4)*2+(sp.spacey/sp.height)];
cb = colorbar;
cb.Position = cbPos;
ax(4).Position = sp.pos{2,2};
ylabel(cb,'PSD (dB re. 1µPa^2/Hz)','interpreter','tex')    
    
    
print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetect4Pan_' moreInfo '.png'],'-dpng','-r150')    
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


% Les isobathes à contourer
isobath = [0:10:150];

% Start mapping
m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');
%m_ruler;

hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');
m_gshhs_f('patch',[.7 .7 .7]); % coastlines

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

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

print([getDirectory('fout') '/generalFigures/presYvanJan/' 'biDetectMap_' moreInfo '.png'],'-dpng','-r150') 
    
end