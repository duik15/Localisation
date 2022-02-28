% This script show the 5 panel diagnostic for a double localisation
% Blue CLD_0950: 19P, 1:56L, 2:28P, 2:44L, 4:41L
% Blue AAV_0950: 16P, 1:56L*, 2:26P, 4:39P
clear all
close all

moreInfo = 'blueSerie';
arrID = {'CLD','AAV'};  
%dtime = [ datetime(2021,07,15, 09, 52, 27); datetime(2021,07,15, 09, 52, 25)];
%dtime = [ datetime(2021,07,16, 06, 19, 03); datetime(2021,07,16, 06, 18, 56)];

dList = [datetime(2021,07,13, 23, 45, 16) datetime(2021,07,13, 23, 45, 13) ;
               datetime(2021,07,13, 23, 45, 35) datetime(2021,07,13, 23, 45, 31) ;
               datetime(2021,07,13, 23, 45, 49) datetime(2021,07,13, 23, 45, 45) ;
               datetime(2021,07,13, 23, 46, 04) datetime(2021,07,13, 23, 45, 59) ;
               datetime(2021,07,13, 23, 46, 22) datetime(2021,07,13, 23, 46, 17) ;
               datetime(2021,07,13, 23, 46, 37) datetime(2021,07,13, 23, 46, 32) ;
               datetime(2021,07,13, 23, 46, 57) datetime(2021,07,13, 23, 46, 54) ;
               datetime(2021,07,13, 23, 47, 09) datetime(2021,07,13, 23, 47, 05) ;
               datetime(2021,07,13, 23, 47, 37) datetime(2021,07,13, 23, 47, 32) ;
               datetime(2021,07,13, 23, 48, 06) datetime(2021,07,13, 23, 48, 02) ;
               datetime(2021,07,13, 23, 48, 23) datetime(2021,07,13, 23, 48, 19) ;
               datetime(2021,07,13, 23, 48, 40) datetime(2021,07,13, 23, 48, 34) ;
               datetime(2021,07,13, 23, 48, 55) datetime(2021,07,13, 23, 48, 50) ;
               datetime(2021,07,13, 23, 49, 34) datetime(2021,07,13, 23, 49, 29) ;
               datetime(2021,07,13, 23, 49, 56) datetime(2021,07,13, 23, 49, 52) ];

azi2 = [57 55 58 57 54 58 57 61 58 58 60 57 60 58 58] - ones(1,size(dList,1))*[2];
azi1 = [69 67 58 71 71 66 66 63 66 66 72 70 63 61 73]  + ones(1,size(dList,1))*[2]; 
aziDBlue = [azi1; azi2]';
           
% Parameter
showFig  = [2 3];
imDur = 15;
buffer = imDur/2;
nbD =1;


% Path 2 folder
%folderWav = ['~/Documents/MPO/BRing/Data/wav/']; 

folderWav = {'/Volumes/BringDD4/Bring_Dep_1/','/Volumes/BringDD3/Bring_Dep_1/'}; 

%folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/honguedo/'];
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/honguedo/' moreInfo '/'];
mkdir2(folderOut)

% Parameter
% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [20 100];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 80];           % [dB] C limite pcolorspgm.im.clims = [30 60];           % [dB] C limite pcolor
%spgm.hyd.clims = [30 80];           % [dB] C limite pcolor
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

d2plot = 4;

for ii=1:2    
dtime(ii) = dList(d2plot,ii)    
%folderWav = folderWav2{ii};    
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

%%
% Calculate PSD and beamForming
for  ii=1:2
ii
folderIn = [folderWav{ii} arrID{ii} '/'];

% Calculate PSD
[~,freqP{ii},timeP{ii},tmp] = spectrogram(wav(ii).pa(:,1),spgm(ii).win.val,spgm(ii).win.novlp,spgm(ii).win.nfft,spgm(ii).fs);
PSD{ii} = 10*log10(tmp);

% Get angle of arrival
%[aziD(ii) detect(ii) ]= fmainBring(arrID{ii},dtime(ii),folderIn);

% Calculate the beamForming
[PSDBF{ii}, timeBF{ii}, freqBF{ii}] = beamForming(arrID{ii}, wav(ii).pa , aziDBlue(d2plot,ii), spgm(ii),'specmethod','spectro');

%if ii==1
%    aziD(1)=240;
%end

end

%%
locD = nan(size(aziDBlue));
for ii=1:size(aziDBlue,1)
    ii
    % And position of location
    locD(ii,:) = beamCrossing(arrID{1},arrID{2}, aziDBlue(ii,1), aziDBlue(ii,2) ,'showFig', 0);
    
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
ylim(spgm(sx).im.freqlims)
set(gca,'FontSize',12)
colormap jet

caxis(spgm(sx).im.clims)




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
ylim(spgm(sx).im.freqlims)
set(gca,'FontSize',12)
colormap jet
caxis(spgm(sx).im.clims)

% Label
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    ylabel('f (Hz)')    
end


xlabel('t (s)')
title(['Beam forming at ' num2str(aziDBlue(sx)) '°'])

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
%locD = beamCrossing(arrID{1},arrID{2}, aziD(1), aziD(2) ,'showFig', 0);

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

print([folderOut 'biDetect5Pan_' moreInfo '.png'],'-dpng','-r150')
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
ylim(spgm(sx).im.freqlims)
set(gca,'FontSize',12)
colormap jet
 
caxis(spgm(sx).im.clims)


%caxis(spgm(ii).im.clims)

% Label
ax(kk).XTickLabel = [];
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    %ylabel('$f$ (Hz)','interpreter','latex')
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
ylim(spgm(sx).im.freqlims)
set(gca,'FontSize',12)
colormap jet
%caxis(spgm(ii).im.clims)
caxis(spgm(sx).im.clims)

% Label
if sx==2; 
    ax(kk).YTickLabel = []; 
else
    %ylabel('$f$ (Hz)','interpreter','latex')    
    ylabel('f (Hz)')    
end


%xlabel('$t$ (s)','interpreter','latex')
xlabel('t (s)')
title(['Beam forming at ' num2str(aziDBlue(d2plot,sx)) '°'])

% Panel index
kk=kk+1;
end

% Colorbar
cbPos = [sp.pos{2,2}(1)+sp.pos{2,2}(3)+0.01 sp.pos{2,2}(2) 0.015 sp.pos{2,2}(4)*2+(sp.spacey/sp.height)];
cb = colorbar;
cb.Position = cbPos;
ax(4).Position = sp.pos{2,2};
ylabel(cb,'PSD (dB re. 1µPa^2/Hz)','interpreter','tex')    
%ylabel(cb,'PSD (dB re. 1$\mu$ Pa${}^2$/Hz)','interpreter','latex')    
    
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
%[locD, distD] = beamCrossing(arrID{1},arrID{2}, aziD(1), aziD(2) ,'showFig', 0);

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
%for ii=1:size(locD,1)
    ii = 4
    colR = ((size(locD,1)-(ii/1.5))/size(locD,1))  *([189, 37, 37]/ 255);
    if ~isnan(locD(ii,1))
        m_plot(locD(ii,2),locD(ii,1),'o','color',colR,'MarkerFaceColor',colR,'MarkerSize',7)
    end
%end

% Ruler
m_ruler([.05 .3],.1,3,'fontsize',8)
%m_northarrow(-65,49.5,1,'type',4,'linewi',.5);

print([folderOut 'biDetectMap_' moreInfo '.png'],'-dpng','-r150') 
    
end