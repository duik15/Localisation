clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file_wav = 'C:\CEDRIC\CHORUS Dropbox\SENSEA_FR\PROJETS\EN_COURS\PROJET_0024_POC_accompagnement _equipe_FV\03_Technique\programme_matlab\AAV_2591LF_20210714T221500_0943_0.wav';
%file_wav = getWavName('')
arrID = 'MLB';
printFig =true;
showFig = [1 2 3 ];
addpath(genpath('../cedBeamForming/'))   % Cedric tools
outName = 'cedVersion_Ns2E14_Voc3_fm150-200hz';
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

[ploc, ptime ]  = getPingLoc(arrID); % Load ping information
%ptime = ptime + seconds(7);

for iFile = 1:length(ptime)
disp(['Loop on file ' num2str(iFile) ' / ' num2str(length(ptime))])
    
close all
[fileName wavID] = getWavName(ptime(iFile), folderIn);
file_wav = [folderIn fileName{1}];
time = ptime(iFile);


% visualisation spectrogramme
%LFFT_spectro = 256;
LFFT_spectro = 2048;
REC = 0.9;
%LFFT_spectro = 2048*8; Baleine bleue
w_pond = kaiser(LFFT_spectro ,0.1102*(180-8.7)); w_pond = w_pond*sqrt(LFFT_spectro/sum(w_pond.^2));
fact_zp =4;

% snapshot sur lequel on va applique la formation de voies "large bande"
[filepath,name] = fileparts(fileName{1});
ftime = getFileTime(name);
start = time - ftime; % Start and end time

t0 = seconds(start) + 2 ; % temps du centre
Ns = 2^16; %65536; % nombre de samples dans le snapshot doit être une puissance de 2

LFFT_FV = Ns;

fmin_int = 150;
fmax_int = 200;

c = 1500;

Lmin = 30;
Lmax = 70;

NB_voies_visees = 72;

SH =-194;
G = 40;
D = 1;

duree_film = 60;
video = ['video_spectro_FV'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
Nc = 20;
R = 10;
theta = 0: 2*pi/(Nc) : 2*pi-2*pi/(Nc);
xc = R*cos(theta);
yc = R*sin(theta);

% chargement des données et dé-normalisation
[MAT_s fe] = audioread(file_wav);
%[MAT_s fe] = wavread(file_wav,[1 1]);
indi_0 = floor(t0*fe);
indi_deb = indi_0-Ns/2+1;
indi_fin = indi_0+Ns/2;

%[MAT_s fe] = wavread(file_wav,[indi_deb indi_fin]);
%[MAT_s fe] = audioread(file_wav,[indi_deb indi_fin]);
file_wav = fileName{1};
% Reading wav
duraNs = Ns / 10000;
[MAT_s,fe, tstart, dura] = readBring([folderIn fileName{1}], ptime(iFile),'duration',duraNs,'buffer',0.4);

MAT_s_vs_t_h = 10^(-SH/20)*10^(-G/20)*D*MAT_s;
aa = size(MAT_s_vs_t_h);
Nsample = aa(1,1);
Nc =aa(1,2);
Duree = (Nsample-1)*1/fe;
t = (0:1:Nsample-1)*1/fe;


% visualisation spectrogramme channel 1
[vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(MAT_s_vs_t_h(:,1),0, fe, LFFT_spectro, REC, w_pond, fact_zp);

figure(1)
pcolor(vec_temps, vec_freq, MAT_t_f_STFT_dB'); shading flat; caxis([30 90]); colorbar
ylim([20 300])
xlabel(' t (s)')
ylabel(' f (Hz)')
title(' Channel 1, dB re. 1µPa2/Hz')
set(gca,'FontSize',16)
colormap jet
caxis([Lmin Lmax])


 % focalisation sur 1 snapshot particulier

 MAT_s_vs_t_h_INT = MAT_s_vs_t_h;

 % formation des FFT des channels
t1 = now;
 for u = 1 : Nc
    MAT_ffts_vs_f_h_INT(:,u)=fft(MAT_s_vs_t_h_INT(:,u));
 end
 vec_f = (0:1:LFFT_FV-1)*fe/LFFT_FV;
 t2 = now;
 temps_formation_FFT_hydro = (t2-t1)*24*3600;
 

 indi_f = find( (vec_f >= fmin_int)&(vec_f <= fmax_int));
 MAT_ffts_vs_f_h_INT_INT = MAT_ffts_vs_f_h_INT(indi_f,:);
 vec_f_INT = vec_f(indi_f);
 

 NN = 360;
 vec_azimut = linspace(0,2*pi-2*pi/NN,NN);
% on raisonne en azimut vs le barycentre du réseau
 x0 = mean(xc);
 y0 = mean(yc);
 xc_rel = xc - x0;
 yc_rel = yc -y0;
% 

t1 = now;
% formation des pondérations frequence, azimut
% cette matrice de pondération pourrait etre calculee en avance et chargé&e
 MAT_POND_vs_azim_h_freq = zeros(length(vec_azimut),length(xc_rel),length(vec_f_INT));
 Ncapt = length(xc_rel);
 for u = 1 : length(vec_azimut)
    for v = 1 : length(xc_rel)
     d_marche =-xc_rel(v)*cos(vec_azimut(u))+-yc_rel(v)*(sin(vec_azimut(u)));
     MAT_POND_vs_azim_h_freq(u,v,:) = 1/Ncapt*exp(-sqrt(-1)*2*pi*vec_f_INT*d_marche/c);
     end
 end
 t2 = now;
 temps_formation_ponderation = (t2-t1)*24*3600;

 
t1 = now;
% goniométrie
for u = 1 : length(vec_azimut)
    MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(u,:,:));
   Energie(u) = 0;
   for v = 1 : length(vec_f_INT)
       vec_pond = (MAT_POND_vs_h_freq(:,v)');
       vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
       Energie(u) = Energie(u)+ (abs(sum(vec_pond.*vec_sfft)))^2;
   end
end
t2 = now;
temps_formation_goniometrie_WB = (t2-t1)*24*3600;


% 
% azimut_NARW = atan2(y_NARW-y0,x_NARW-x0);
% azimut_SHIP = atan2(y_SHIP-y0,x_SHIP-x0);
Energie_dB = 10*log10(Energie);
figure(2)
plot(vec_azimut*180/pi, Energie_dB,'k','LineWidth',2);
hold on
xlabel(' azimut (°)')
ylabel(' SPL dB re. 1µPa')
grid on
set(gca,'FontSize',16)
% for u = 1 : length(azimut_NARW)
%    [val indi] = min(abs( vec_azimut - azimut_NARW(u)));
%    plot(vec_azimut(indi)*180/pi,Energie_dB(indi),'bs','MarkerFaceColor','b')
% end
% for u = 1 : length(azimut_SHIP)
%    [val indi] = min(abs( vec_azimut - azimut_SHIP(u)));
%    plot(vec_azimut(indi)*180/pi,Energie_dB(indi),'r^','MarkerFaceColor','r')
% end


figure(3)
R_source = 5000;
vec_azim = linspace(0,2*pi,1000);
x_cercle = R_source*cos(vec_azim);
y_cercle = R_source*sin(vec_azim);
plot(x0,y0,'k+')
hold on
plot(xc,yc,'r.','MarkerSize',12)
plot(x_cercle,y_cercle,'k')
% plot(x_NARW,y_NARW,'^b','MarkerFaceColor','b','MarkerSize',14);
% plot(x_SHIP,y_SHIP,'sk','MarkerFaceColor','k','MarkerSize',14);
xlabel(' x (m)')
ylabel(' y (m)')
grid on
axis equal
for u = 1 : length(xc)
   a = atan2(yc(u),xc(u));
   xx = R_source*cos(a);
   yy = R_source*sin(a);
   plot(xx,yy,'r.');
   text(xx,yy,['h' num2str(u)])
end
% for u = 1 : length(x_NARW)
%    text(x_NARW(u),y_NARW(u),['NARW' num2str(u)])
% end
% for u = 1 : length(x_SHIP)
%    text(x_SHIP(u),y_SHIP(u),['SHIP' num2str(u)])
% end

Energie_NORM = Energie/max(Energie);
xxx = cos(vec_azimut).*Energie_NORM*R_source;
yyy = sin(vec_azimut).*Energie_NORM*R_source;
plot(xxx,yyy,'k','LineWidth',2)

 % formation d'une voie 
azimut_cible =30;
[ val indi_azimut_vise] = min(abs(vec_azimut - azimut_cible));
t1 = now;
 MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(indi_azimut_vise,:,:));
   for v = 1 : length(vec_f_INT)
       vec_pond = (MAT_POND_vs_h_freq(:,v)');
       vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
       vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
   end
FFT_S_FV_int = zeros(1,LFFT_FV/2+1);
FFT_S_FV_int(indi_f) = vec_sfft_FV;
FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
s_FV = ifft(FFT_S_FV);
t2 = now;
temps_formation_1voie = (t2-t1)*24*3600;
% spectrogramme 1 voie
figure(7)
t = t0-Ns/2*1/fe+(0:Ns-1)*1/fe;
subplot(5,1,1);plot(t,s_FV,'k')
xlabel(' t (s)');
[vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
subplot(5,1,[2:5])
pcolor(vec_temps_FV, vec_freq_FV,MAT_t_f_STFT_dB_FV'); shading flat; caxis([30 90]); colorbar
ylim([fmin_int fmax_int])
xlabel(' t (s)')
ylabel(' f (Hz)')
title(['VOIE n°' num2str(indi_azimut_vise) ', dB rel'])
caxis([Lmin Lmax])
colormap jet


% formation de N voies
t1 = now;
 spie = 1;
 pas = floor(length(vec_azimut)/NB_voies_visees);
 MAT_s_FV_vs_t_voie=zeros(Ns,length(1 : pas : length(vec_azimut)));
 
 for u = 1 : pas : length(vec_azimut)
     %u
     vec_azim_voie(spie) = vec_azimut(u);
     indice_voie(spie) = u;
%     
% 
 MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(u,:,:));
    for v = 1 : length(vec_f_INT)
        vec_pond = (MAT_POND_vs_h_freq(:,v)');
        vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
        vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
    end
 FFT_S_FV_int = zeros(1,LFFT_FV/2+1);
 FFT_S_FV_int(indi_f) = vec_sfft_FV;
 FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
 s_FV = ifft(FFT_S_FV);
% 
 MAT_s_FV_vs_t_voie(:,spie)=s_FV;
 spie = spie + 1;
 end
t2 = now;
temps_formation_voies_WB = (t2-t1)*24*3600;

% calcul des spectrogrammes des voies
t1 = now;
Nvoie = length(vec_azim_voie);
vec_freq_FV1 =(0:1:LFFT_spectro*fact_zp-1)*fe/(LFFT_spectro*fact_zp);
indi_f = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
 MAT_S_FV_vs_voie_t_f = zeros(Nvoie,length(vec_temps_FV),length(indi_f));
 for u = 1 : Nvoie
%     u
%     Nvoie
     [vec_temps_FV, vec_freq_FV1, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(MAT_s_FV_vs_t_voie(:,u),t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
  indi_f = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
  MAT_S_FV_vs_voie_t_f(u,:,:)=MAT_t_f_STFT_dB_FV(:,indi_f);
%    
end
 vec_freq_FV = vec_freq_FV1(indi_f);
 t2 = now;
 temps_formation_spectro_voie = (t2-t1)*24*3600;

 % tracé des spectrogrammes
 %{
 t1 = now;
 figure(8)
 a =sqrt(Nvoie)
 if a == floor(a)
     a=a;
 else
     a=floor(a)+1;
 end
 
 
writerObj = VideoWriter(video) 
writerObj.FrameRate = (Nvoie)/duree_film ;
writerObj.Quality = 100;
open(writerObj)

for u = 1 : Nvoie
%      u
%      Nvoie
figure(4)
    M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
    pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat; caxis([30 90])
    ylim([fmin_int fmax_int])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
 title(['V ' num2str(u) 'Az (°) ' num2str(floor(180/pi*vec_azim_voie(u)))])
 caxis([Lmin Lmax])
 colorbar
 colormap jet
 
     
h = gcf;
F = getframe(h);
writeVideo(writerObj,F)
close
end

  close(writerObj)
t2 = now;

 temps_affichage_spectro_voie = (t2-t1)*24*3600;

 
 figure(5)
 for u = 1 : 1: Nvoie/2
%      u
%      Nvoie
    subplot(6,6,u)
    M =  squeeze(MAT_S_FV_vs_voie_t_f(2*u,:,:));
    pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat;
    ylim([fmin_int fmax_int])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
 title(['V ' num2str(2*u) 'Az (°) ' num2str(floor(180/pi*vec_azim_voie(2*u)))])
 caxis([Lmin Lmax])
 colormap jet
 end
% 
% 
%}

if printFig ==true
   figName = {'spectroChanel1', 'energyFctAzimutXY','energyFctAzimutCircle',...
              'spectroBestAngle','spectroNAngle', 'subplotSpectroNAngle'};
    for jj=1:length(showFig)
       fid = ['-f' num2str(showFig(jj))];
       print([folderOut figName{showFig(jj)} '_' outName '_p' num2str(iFile)  '_fID' wavID{1} '.png'], '-r150','-dpng', fid) 
   end
    
end

end  % End loop