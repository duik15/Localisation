% This script run beamforming over some ping emission to find back the
% array.
% This version as the sin(X) - cos(y) for the xc,yc but not for all other cos,sin.
clear all
close all


% Path information
arrID = 'PRC';
addpath(genpath('../cedBeamForming/'))   % Cedric tools
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/'];
%folderIn = 'Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\';
outName = 'test_xysc_OffSet135_voc3_Ns2e14_prcCircle';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];

% Loading files and time
% Or give an array of time to process ptime = [datetime(xxx),datetime(xxx)]
[ploc, ptime ]  = getPingLoc(arrID); % Load ping information
ptime = ptime + seconds(7);
%ptime =  datetime(2021,07,15,14,32,50);
[fileList wavID] = getWavName(ptime, folderIn);
wavInfo = audioinfo([folderIn fileList{1}]);
nbF = length(fileList);

% Figure parameters
saveData = false;
showFig = [3 7]%1 2 3 4 5 6 7]; %7 9]%[1 2 3 4 7];% [1:7];% Figure number to print
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^14;              % Total number of sample
buffer = 0.4;             % Time in second
duraNs = Ns / 10000;      % Length of wav in second to load

if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Spectro parameter
% visualisation spectrogramme

% Spectro windows
LFFT_spectro = 2048;
REC = 0.9;
w_pond = kaiser(LFFT_spectro ,0.1102*(180-8.7)); w_pond = w_pond*sqrt(LFFT_spectro/sum(w_pond.^2));
fact_zp =4;

% Frequence min and max
fmin_int = 150;
fmax_int = 200;

% Sound velocities
c = 1500;

% Pcolor plot
Lmin = 30;
Lmax = 70;

NB_voies_visees = 72;

SH =-194;
G = 40;
D = 1;

duree_film = 60;

% Creating a circle vecteur
Nc = 20;
R = 10;
theta = 0: 2*pi/(Nc) : 2*pi-2*pi/(Nc);
offset = 135 * pi /180; %(90 + 25) * pi / 180;  % Av = (90-67) CD  =
xc = R*sin(theta +  offset);
yc = R*cos(theta +  offset);



%% File loop

% Initialisation
matEnergie = nan(nbF,360);
angleA = nan(nbF,nbPk);

for iFile =1:length(fileList)
    close all
    
    %iFile = 3
    disp(['Executing beamforming for file ' num2str(iFile) '/' num2str(length(fileList))])
    % Loading part
    file_wav = fileList{iFile};
    % Reading wav
    [MAT_s,fe, tstart, dura] = readBring([folderIn fileList{iFile}], ptime(iFile),'duration',duraNs,'buffer',buffer,'power2',true);
    
    %Ns = dura * fe;
    LFFT_FV = Ns;
    
    
    MAT_s_vs_t_h = 10^(-SH/20)*10^(-G/20)*D*MAT_s;
    aa = size(MAT_s_vs_t_h);
    Nsample = aa(1,1);
    Nc =aa(1,2);
    Duree = (Nsample-1)*1/fe;
    t = (0:1:Nsample-1)*1/fe;
    
    
    % visualisation spectrogramme channel 1
    [vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(MAT_s_vs_t_h(:,1),0, fe, LFFT_spectro, REC, w_pond, fact_zp);
    
    
    % focalisation sur 1 snapshot particulier
    
    MAT_s_vs_t_h_INT = MAT_s_vs_t_h;
    
    % formation des FFT des channels
    t1 = now;
    MAT_ffts_vs_f_h_INT = nan(size(MAT_s_vs_t_h_INT));
    
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
            d_marche =-xc_rel(v)*sin(vec_azimut(u))+-yc_rel(v)*(cos(vec_azimut(u)));
            MAT_POND_vs_azim_h_freq(u,v,:) = 1/Ncapt*exp(-sqrt(-1)*2*pi*vec_f_INT*d_marche/c);
        end
    end
    t2 = now;
    temps_formation_ponderation = (t2-t1)*24*3600;
    
    
    % Goniométrie -------------------> Figure 2,3
    t1 = now;
    Energie = [];
    for u = 1 : length(vec_azimut)
        MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(u,:,:));
        Energie(u) = 0;
        for v = 1 : length(vec_f_INT)
            vec_pond = (MAT_POND_vs_h_freq(:,v)');
            vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
            Energie(u) = Energie(u)+ (abs(sum(vec_pond.*vec_sfft)))^2;
        end
    end
    
    % Other form of Energie unity
    Energie_NORM = Energie/max(Energie);
    Energie_dB = 10*log10(Energie);
    
    % Find peak of energy
    [pk pkloc] = findpeaks(Energie_NORM,'SortStr','descend','NPeaks', 4 );%'MinPeakHeight',0.3);
    
    % Find the direction of source
    angleA(iFile, : )  = pkloc;  % Angle of arrival with side lobe
    angleM(1,iFile)  = pkloc(1);  % Max angle of arrival

    

    % Saving output multi figure
    matEnergie(iFile, :) = Energie;
    
    
    
 % ---------------------------------------------------------------------   
 % ------------------------ Figure in loop -----------------------------    
    % -------------------------------------------------------------------
% Figure 1 : spectogram
if any( showFig == 1 )
    figure(1)
    pcolor(vec_temps, vec_freq, MAT_t_f_STFT_dB'); shading flat; caxis([30 90]); colorbar
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1µPa2/Hz')
    set(gca,'FontSize',16)
    colormap jet
    %caxis([Lmin Lmax])
end


% -------------------------------------------------------------------
% Figure 2 : Energy fct azimut
if any( showFig == 2 )
    % azimut_NARW = atan2(y_NARW-y0,x_NARW-x0);
    % azimut_SHIP = atan2(y_SHIP-y0,x_SHIP-x0);
    
    figure(2)
    plot(vec_azimut*180/pi, Energie_dB,'k','LineWidth',2);
    hold on
    xlabel(' azimut (°)')
    ylabel(' SPL dB re. 1µPa')
    grid on
    set(gca,'FontSize',16)
end

% -----------------------------------------------------------------
% Figure 3 : Energy fct azimut circle
if any( showFig == 3 )
    figure(3)
    R_source = 5000;
    vec_azim = linspace(0,2*pi,1000);
    x_cercle = R_source*sin(vec_azim);
    y_cercle = R_source*cos(vec_azim);
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
        a = atan2(xc(u),yc(u));
        xx = R_source*sin(a);
        yy = R_source*cos(a);
        plot(xx,yy,'r.');
        text(xx,yy,['h' num2str(u)])
    end
    % for u = 1 : length(x_NARW)
    %    text(x_NARW(u),y_NARW(u),['NARW' num2str(u)])
    % end
    % for u = 1 : length(x_SHIP)
    %    text(x_SHIP(u),y_SHIP(u),['SHIP' num2str(u)])
    % end
    
   
    xxx = sin(vec_azimut).*Energie_NORM*R_source;
    yyy = cos(vec_azimut).*Energie_NORM*R_source;
    plot(xxx,yyy,'k','LineWidth',2)
end % end showfig3

% ---------------------------------------------------------------------
% Figure 4: spectrogramme 1 voie
if any ( showFig == 4 )
    
    % formation d'une voie
    [~, iamax ]  =max(Energie_dB);
    azimut_cible = vec_azimut(iamax);
    [ val indi_azimut_vise] = min(abs(vec_azimut - azimut_cible));
    t1 = now;
    MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(indi_azimut_vise,:,:));
    for v = 1 : length(vec_f_INT)
        vec_pond = (MAT_POND_vs_h_freq(:,v)');
        vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
        vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
    end
    FFT_S_FV_int = zeros(1,floor(LFFT_FV/2)+1);
    FFT_S_FV_int(indi_f) = vec_sfft_FV;
    FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
    s_FV = ifft(FFT_S_FV);
    t2 = now;
    temps_formation_1voie = (t2-t1)*24*3600;
     
    figure(4)
    t0=1;
    t =t0-Ns/2*1/fe+(0:Ns-1)*1/fe;
    subplot(5,1,1);plot(t,s_FV,'k')
    xlabel(' t (s)');
    %[vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(MAT_s_vs_t_h(:,1),0, fe, LFFT_spectro, REC, w_pond, fact_zp);
    [vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
    subplot(5,1,[2:5])
    pcolor(vec_temps_FV, vec_freq_FV,MAT_t_f_STFT_dB_FV'); shading flat; caxis([30 90]); colorbar
    ylim([fmin_int fmax_int])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(['VOIE n°' num2str(indi_azimut_vise) ', dB rel'])
    %caxis([Lmin Lmax])
    colormap jet
end

% ---------------------------------------------------------------
% Figure 5 : spectogram de N voies
if any (showFig ==5 )
    
    % 
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
        FFT_S_FV_int = zeros(1,floor(LFFT_FV/2)+1);
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
    t1 = now;
   
    a =sqrt(Nvoie);
    if a == floor(a)
        a=a;
    else
        a=floor(a)+1;
    end
    
    videoName = ['spectro_' outName '_p' num2str(iFile)  '_fID' wavID{iFile}  ];
    writerObj = VideoWriter([folderOut videoName '.mp4'],'MPEG-4');
    writerObj.FrameRate = (Nvoie)/duree_film ;
    writerObj.Quality = 100;
    open(writerObj)
    
    for u = 1 : Nvoie
        %      u
        %      Nvoie
        figure(5)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat; caxis([30 90])
        ylim([fmin_int fmax_int])
        xlabel(' t (s)')
        ylabel(' f (Hz)')
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(180/pi*vec_azim_voie(u)))])
        caxis([Lmin Lmax])
        colorbar
        colormap jet
        
        %print([folderOut 'spectroNAngle_' outName '_' num2str(vec_azim_voie(u)* 360 / 2/pi)   'd.png'], '-r150','-dpng', '-f5') 
        
        h = gcf;
        F = getframe(h);
        writeVideo(writerObj,F)
        close
        
    end
    
    close(writerObj)
    t2 = now;
    
    temps_affichage_spectro_voie = (t2-t1)*24*3600;
end % end showfig 4


%  -------------------------------------------------------------------
% Figure 6 :  Spectogram N vois 1 figure
if any ( showFig == 6)
    figure(6)
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
    
    
end

% Saving figure
if printFig ==true
   figName = {'spectroChanel1', 'energyFctAzimutXY','energyFctAzimutCircle',...
              'spectroBestAngle','spectroNAngle', 'subplotSpectroNAngle'};
    for jj=1:length(showFig)
       fid = ['-f' num2str(showFig(jj))];
       if showFig(jj) ~= 5 && showFig(jj) ~= 7 && showFig(jj)<7
       print([folderOut figName{showFig(jj)} '_' outName '_p' num2str(iFile)  '_fID' wavID{iFile} '.png'], '-r150','-dpng', fid) 
       end
   end
    
end

% -----------------------------------------------------------------------    

end % end loop on file



%% Saving data
if saveData == true
save([folderOut 'dataAngleOfArrival_' outName '.mat'],'angleA','angleM','ptime','ploc','matEnergie')
end
%% ---------------- Figures -------------------------------


if any ( showFig == 7)

if sum(diff(angleM(3:end-2)))< 0
    pNb = 10:-1:1;
    angleA = fliplr(aMax3);
else
    pNb = 1:10;
end    
    
% Linéaire regression
x = 1:10;
p = polyfit(x,angleM,1);
angleFit = polyval(p,x);  


figure(7)
hA = plot(x,angleM,'-.ok');
hold on
hTh = plot(x,0:36:360-36,'k-');
hFit = plot(x,angleFit);

% Plot some side lobe
for ii=2:2
   hSl = plot(x,angleA(:,ii),'o') 
end

xlabel('Ping number')
ylabel('Angle of arrival')
ax = gca;
ax.XTickLabel = pNb;
leg= legend([hA,hTh,hFit, hSl],{'Ping','Theoritical curve','Lin. Reg. Curve','Secondary lobe'},'Location','northwest')

print([folderOut 'angleOfArrival_' outName '.png'],'-r150','-dpng')

end

if any ( showFig ==8)
figure(8)

plot(xc,yc,'o')


for ic = 1:length(xc)
text(xc(ic)+0.2,yc(ic),num2str(ic),'fontsize',14)
end
end

if any ( showFig == 9)
    
    
    figure(9)
    
    plot(aMax,'--ok')
    xlabel('Ping number')
    ylabel('Angle of arrival')
    print([folderOut 'angleOfArrivalbad_' outName '.png'],'-r150','-dpng')
    
    
    
     [asort, isort] = sort(aMax);
    
    figure(7)
    
    plot(asort,'--ok')
    hold on
    plot([])
    xlabel('Ping number')
    ylabel('Angle of arrival')
    
    % Thick sort
    ax=gca;
    
    vec10 = 1:10;
    pingSort = vec10(isort)
    pingSortStr = strsplit(num2str(vec10(isort)));
    ax.XTickLabel = pingSortStr;
    print([folderOut 'angleOfArrival_sort_' outName '.png'],'-r150','-dpng')
end