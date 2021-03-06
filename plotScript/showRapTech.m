% This script show figure for technical repport
close all

% figure 1 : Show 10 chanel
% figure 2 + 3 + 4: Generate 6 differents angles + gonio clean  + map
% figure 5 : ?? Localisation
% figure 6 + 7 : Theorique circle + theorique signal

%% ---------------------------------------------------------------

% figure 1 : Show 10 chanel
sp.height=30; sp.width=30; sp.nbx=2; sp.nby=5;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.2; sp.spacey=0.2;
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');
%mesh2= mesh2';
nbSP = sp.nbx * sp.nby;
fontS = 12;

i_ch = 1:2:2*nbSP;

for ii=1:nbSP
    [~,freq,time,tmp] = spectrogram(wav.pa(:,i_ch(ii)),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    if ii==1; PSD2 = nan([nbSP size(tmp)]) ; end
    PSD2(ii,:,:) = 10*log10(tmp);
end

fig = figure(1);
movegui(fig,'northwest');
% --------------- Figure --------------------------
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [sp.width sp.height]);
set(fig, 'PaperPositionMode', 'manual');

for ii=1:nbSP
    [sx,sy] = find(mesh2'==ii);
    % Get the spectogram of first chanel
    
    PSD = squeeze(PSD2(ii,:,:));
    ax(ii)=axes('position',sp.pos{sx,sy},'FontSize',fontS,'Box','on','Layer','top');
    pcolor(time,freq,PSD);
    shading flat;
    colormap(jet)
    %axis tight;
    ylim(spgm.im.freqlims);
    %set(gca,'YScale',spgm.im.freqscale);
    caxis([spgm.im.clims])

    
    % Text and label
    if ii<=5
        text(1.68, spgm.im.freqlims(2) * 0.93,['ch.' num2str(i_ch(ii))], 'backgroundcolor', [1 1 1],'margin',2)
    else
        text(1.63, spgm.im.freqlims(2) * 0.93,['ch.' num2str(i_ch(ii))], 'backgroundcolor', [1 1 1],'margin',2)
    end
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2 + 5,['ch.'])
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2-5,[num2str(mainCh(ii))])
    if sx == 1
        yl = ylabel('f (Hz)','FontSize',fontS);
    else
        ax(ii).YTickLabel = [];
    end
    if sy == sp.nby
        xl = xlabel('t (s)','FontSize',fontS) ;
    else
        ax(ii).XTickLabel = [];
    end
    
    set(gca,'layer','top','XMinorTick','on','YMinorTick','on');
    grid on
    box on
    
end



% Color bar
cbPos = [sp.pos{sx,sy}(1)+sp.pos{sx,sy}(3)+0.01 sp.pos{sx,sy}(2) 0.02 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
cb.Position(3) = 0.02;
%ax(end).Position = sp.pos{end};
ylabel(cb,'PSD (dB re. 1?Pa^2/Hz)')

% General text and label
%sgtitle([arrID ' | ' datestr(ptime)],'Interpreter','none');

print('-dpng','-r150',[folderOut 'fig1_rapTech_showChanels.png'])



%% ---------------------------------------------------------------
% Figure 2 : Generate 6 differents angles

clear ax
sp.height=10; sp.width=10; sp.nbx=1; sp.nby=1;
sp.ledge=0.2; sp.redge=0.2; sp.tedge=0.2; sp.bedge=0.2;
sp.spacex=0.3; sp.spacey=0.3;
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');
fontS = 12;

nbFig =  9;

% Azimut to plot
azimutV = 0:360/((nbFig)-1):360;

% Get the beamforming spectrogram
[Pdb, timeV, freqV] = beamForming(arrID, wav.pa , azimutV, spgm,'specmethod','spectro');

for ii=1:nbFig
    
    fig = figure(2);
    ax = axes('position',sp.pos{1},'FontSize',fontS,'Box','on','Layer','top');
    set(fig, 'PaperPosition', [0 0 5 5]);    % can be bigger than screen
    M =  squeeze(Pdb(ii,:,:));
    pcolor(timeV, freqV,M'); shading flat;
    ylim([spgm.im.freqlims])
    caxis([30 70])
    %ylabel(' f (Hz)')
    %xlabel(' t (s)')
    
    set(gca,'XTickLabel',[],'YTickLabel',[],'layer','top')
    colormap jet
    grid on
    box on
    axis tight
    print([folderOut 'fig2_rapTech_snapshotVoie_azimut' num2str(azimutV(ii)) '.png'],'-dpng','-r150')
    clf
end
%%
clear ax
sp.height=10; sp.width=10; sp.nbx=1; sp.nby=1;
sp.ledge=1.3; sp.redge=0.2; sp.tedge=0.2; sp.bedge=1.1;
sp.spacex=0.3; sp.spacey=0.3;
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');
fontS = 12;

fig = figure(2);
ax = axes('position',sp.pos{1},'FontSize',fontS,'Box','on','Layer','top');
set(fig, 'PaperPosition', [0 0 5 5]);    % can be bigger than screen
M =  squeeze(Pdb(2,:,:));
pcolor(timeV, freqV,M'); shading flat;
ylim([spgm.im.freqlims])
caxis([30 70])
ylabel(' f (Hz)')
xlabel(' t (s)')
colormap jet
set(gca,'layer','top')
grid on
box on
axis tight
print([folderOut 'fig2_rapTech_snapshotVoie_azimut_tick_' num2str(azimutV(ii)) '.png'],'-dpng','-r150')


%% ---------------------------------------------------------------
% Figure 3 : Gonio clean

aziCircle = linspace(0,2*pi,1000);
arrDiameter = 10;

% Theorique hydro position
Nc = 20;
R = 10;
theta = 0: 2*pi/(Nc) : 2*pi-2*pi/(Nc);
xc = R*sin(theta +  arr.offSet * pi /180);
yc = R*cos(theta +  arr.offSet * pi /180);
xt = (R-1.5)*sin(theta +  arr.offSet * pi /180) -1;
yt = (R-1.5)*cos(theta +  arr.offSet * pi /180);

fig = figure(3);
set(fig, 'PaperPosition', [0 0 10 10]);    % can be bigger than screen

plot(mean(arr.xh),mean(arr.yh),'k+')
hold on
plot(arrDiameter*sin(aziCircle),arrDiameter*cos(aziCircle),'k')
plot(xc,yc,'r.','MarkerSize',18)
for u = 1 : length(arr.xh)
    text(xt(u)+0.3,yt(u),[ num2str(u)],'fontsize',16)
end

% Plot energie norm
xxx = sin(azimut360* pi /180).*Energie_NORM*arrDiameter;
yyy = cos(azimut360* pi /180).*Energie_NORM*arrDiameter;
plot(xxx,yyy,'k','LineWidth',2)

% Label, axes and title
xlabel(' x (m)')
ylabel(' y (m)')
xlim( [-15 15]);
ylim( [-15 15]);
axis equal

% Draw a north arrow
%quiver(9,7,0,3,'Color',[0 0 0],'LineWidth',1)
%text(9.2,9.5,'N')


print([folderOut 'fig2_rapTech_goniogram.png'],'-dpng','-r150')

%% ---------------------------------------------------------------
% Figure 4 : Map
clf
% Load map
gsl = load([getDirectory('map') 'GSL_bathy_500m.mat']);
gebco = load([getDirectory('map') 'gebco_colormap.dat']);

[mlbLoc, mlb]= getArrInfo('MLB');
[prcLoc, prc]= getArrInfo('PRC');

% Default parameter
winSize = 0.05;

% Zone de la carte
lon_min = -64.4;
lon_max = -63.9;
lat_min =  48.4;
lat_max =  48.7;

% Les isobathes ? contourer
isobath = [0:10:450];

fig = figure(4);
set(fig, 'PaperPosition', [0 0 7 7]);    % can be bigger than screen

m_proj('mercator', 'long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_grid('fontsize',12,'linestyle','none');

% Un beau colormap pour la bathymetrie qui vient de GEBCO
colormap(gebco);

hold on
m_contourf(gsl.LON,gsl.LAT,gsl.Z,isobath,'linestyle','none');

% Information graphique
c = colorbar('location','eastoutside');
cpos = c.Position;
set(c, 'FontSize', 12, 'FontName','times')
ylabel(c,'Depth(m)')

xlabel('Longitude', 'FontSize', 12, 'FontName','times')
ylabel('Latitude', 'FontSize', 12, 'FontName','times')

m_gshhs_f('patch',[.7 .7 .7]); % coastlines


% ------------ Plot points ----------------
hm=   m_plot(mlbLoc(2),mlbLoc(1),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor','r','MarkerSize',6) ;
hm=   m_plot(prcLoc(2),prcLoc(1),'o','MarkerEdgeColor','r',...
    'MarkerFaceColor','r','MarkerSize',6) ;

print([folderOut 'fig2_rapTech_map.png'],'-dpng','-r150')


%% -------------------------------------------------------------
% figure 7 : Theorique signal
% Parameter
nbZ = 250;
yOff = 0;
% Sin variable
tt = linspace(0,8,100);
nbT = 2* nbZ + 100;
%tt2 = linspace(0,16,200);
omega = [ 2 3 4 ];
amp = [ 10 1 1 ];
istart = floor(nbT/2 - nbZ *  [1/4  1/2 3/4 ]  ) ;
noiseAmp = 0.3;

% Get color order
colorO = get(gca,'colorOrder');
clf

% Create the sin matrice
sin0= zeros(1,length(tt));
for ii=1:length(omega)
    sin0 = sin0 + amp(ii) * sin(omega(ii) .* tt  + pi/2);
end

% Gaussian function
sig = 1;
gau  = -10*(2/sqrt(2*pi*sig^2)) * exp(-(tt - (tt(end)/2)).^2/(2* sig^2));
sin2 =  [zeros(1,nbZ) sin0 .* gau zeros(1,nbZ)];
sin2 = sin2 / max(sin2);
sin2 = sin2 - yOff;

sinM = zeros(length(istart),nbZ);
for io=1:length(istart)
    sinR(io,:) = sin2 + noiseAmp*( rand(1,length(sin2)) - 1/2);
    sinM(io,:) = sinR(io,istart(io):istart(io)+nbZ-1);
end


%Subplot parameter
clear ax
sp.height=10; sp.width=10; sp.nbx=1; sp.nby=3;
sp.ledge=0.9; sp.redge=0.2; sp.tedge=0.2; sp.bedge=1.3;
sp.spacex=0.3; sp.spacey=0.2;
[sp.pos, mesh2]=subplot2(sp,'meshorder','y');
fontS = 12;

%% ----------------------------------------------------- 
% Fig4 : Not Sync Signal

fig = figure(7);
%colorO = get(gca,'colorOrder');
set(fig, 'PaperPosition', [0 0 7 4]);    % can be bigger than screen


for io=1:length(istart)
    [sx,sy] = find(mesh2'==io);
    ax(io)=axes('position',sp.pos{sx,sy},'FontSize',fontS,'Box','on','Layer','top');
    hold on
    %plot([0 nbZ], [0 0] ,'-.k')
    plot(sinM(io,:) ,'color',colorO(io+1,:))
    plot(floor(nbZ/4 * io * [1 1]) , [-1.5 1.5],'-.','color', colorO(io+1,:))
    hold off

ax=gca;

ax.XTick = [0:nbZ/4:nbZ];
if io==3
xlabel('Time')
ax.XTickLabel = {'0','1t_0','2t_0','3t_0','4t_0'};
else
    ax.XTickLabel = [];
end
if io ==2
ylabel('Amplitude')
end
ax.YTick = [-1:1:1];
ax.YTickLabel = strsplit(num2str([-1:1:1]));
xlim([0 nbZ])
ylim([-1.15 1.15])
box on
grid on

%ax.GridColor = 'k'
ax.GridAlpha = 0.5 % maximum line opacity

end

print([folderOut 'fig4_rapTech_theoSignal_notSync.png'],'-dpng','-r150')

%%

fig = figure(7);
%colorO = get(gca,'colorOrder');
set(fig, 'PaperPosition', [0 0 7 4]);    % can be bigger than screen

for io=1:length(istart)   
    [sx,sy] = find(mesh2'==io);
    ax(io)=axes('position',sp.pos{sx,sy},'FontSize',fontS,'Box','on','Layer','top');
    %plot([0 nbZ], [0 0] ,'-.k')
    plot(sinR(io,175 :175 + nbZ-1) ,'color',colorO(io+1,:))

    ax=gca;
    
    ax.XTick = [0:nbZ/4:nbZ];
    if io==3
        xlabel('Time')
        ax.XTickLabel = {'0','1t_0','2t_0','3t_0','4t_0'};
    else
        ax.XTickLabel = [];
    end
    if io==2
    ylabel('Amplitude')
    end
    ax.YTick = [-1:1:1];
    ax.YTickLabel = strsplit(num2str([-1:1:1]));
    xlim([0 nbZ])
    ylim([-1.15 1.15])
    box on
    grid on
    ax.GridAlpha = 0.5 % maximum line opacity
end

print([folderOut 'fig4_rapTech_theoSignal_Sync.png'],'-dpng','-r150')

%%

% Sum all singal
sumSin = sum(sinR,1);

clf
fig = figure(7);
colorO = get(gca,'colorOrder');
set(fig, 'PaperPosition', [0 0 7 4]);    % can be bigger than screen
hold on
plot([0 nbZ], [0 0] ,'-.k')
plot(sumSin(175 :175 + nbZ-1),'k')
hold off

ax=gca;
xlabel('Time')
ylabel('Amplitude')
ax.XTick = [0:nbZ/4:nbZ];
ax.XTickLabel = {'0','1t_0','2t_0','3t_0','4t_0'};
ax.YTick = [-5:1:5];
ax.YTickLabel = strsplit(num2str([-5:1:5]));
xlim([0 nbZ])
ylim([-3.3 3.3])
box on
grid on

ax.GridAlpha = 0.5 % maximum line opacity

print([folderOut 'fig4_rapTech_theoSignal_sumsyn.png'],'-dpng','-r150')

%% Sum destructive


% Sum all singal

%sumSinD = zeros(size(sinR,2));
%for io=1 : size(sinR,1)
%    sumSinD = sumSinD + sinR
%end
    
clf
fig = figure(7);
colorO = get(gca,'colorOrder');
set(fig, 'PaperPosition', [0 0 10 5]);    % can be bigger than screen
hold on
plot([0 nbZ], [0 0] ,'-.k')
%plot(sum(sinM,1),'k')
plot(sinM,'k')
hold off

ax=gca;
xlabel('Time')
ylabel('Amplitude')
ax.XTick = [0:nbZ/5:nbZ];
ax.XTickLabel = {'0','1t_0','2t_0','3t_0','4t_0','5t_0'};
ax.YTick = [-5:1:5];
ax.YTickLabel = strsplit(num2str([-5:1:5]));
xlim([0 nbZ])
ylim([-5 5])
box on

print([folderOut 'fig4_rapTech_theoSignal_sumsinD.png'],'-dpng','-r150')




%% Figure 8 : Gognio for time

arrDiameter = 10;
offset = 0;
tethaCircle = linspace(0,2*pi,100);
% Theorique hydro position
Nc = 20;
R = 10;
theta = 0: 2*pi/(Nc) : 2*pi-2*pi/(Nc);
xc = R*sin(theta +  offset * pi /180);
yc = R*cos(theta + offset * pi /180);
xt = (R-1.5)*sin(theta +  offset * pi /180) -1;
yt = (R-1.5)*cos(theta +  offset * pi /180);

fig = figure(3);
set(fig, 'PaperPosition', [0 0 10 10]);    % can be bigger than screen
plot(mean(xc),mean(yc),'k.','markersize',12)
hold on
plot(arrDiameter*sin(tethaCircle),arrDiameter*cos(tethaCircle),'k')
plot(xc,yc,'r.','MarkerSize',20)
for u = 1 : length(arr.xh)
    text(xt(u)+0.3,yt(u),[ num2str(u)],'fontsize', 18)
end

% Plot time line
for it=1:4
   t0 = 20 - (10* (it-1))   ;
   plot([-12.2 12], t0 * [1 1],'--','color',colorO(it,:))
   text(-12,t0-1,[num2str(it-1) 't_0'],'color',colorO(it,:),'fontsize', 20)
end
%annotation('textarrow',0.3*[1 1], [0.83 0.2],'String','','FontSize',13,'Linewidth',2)
%quiver(-12,20,0,-34,1,'k')
text(-11.2,-13,'t','fontsize',20)

xlim([-20 20])
ylim([-15 25])
axis square

annotation('textarrow',0.282*[1 1], [0.807 0.17],'String','','FontSize',13,'Linewidth',2)
print([folderOut 'fig4_rapTech_goniogram.png'],'-dpng','-r150')
%%


