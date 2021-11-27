% This is just a script containing the figure for each timestepd on locateBring
% ---------------------------------------------------------------------
% Figure 1 : spectogram
% Figure 2 : Energy fct azimut
% Figure 3 : Energy fct azimut circle
% Figure 4: spectrogramme 1 voie
% Figure 5 :  Spectogram N vois that cover 360deg 1 figure
% Figure 6 :  Spectogram N vois around azi cible 1 figure
% Figure 9 : spectogram de N voies in a forloop for video

% ------------------- General Parameter ------------------------------
% formation de voie
if ~exist('aziCible') || strcmp(aziCible,'max')
    aziCible = angleM;
end
[ val indAziCible] = min(abs(vec_azimut*180/pi - aziCible));



% ------------------------ Figure in loop -----------------------------

% Figure 1 : spectogram
if any( showFig == 1 )
    figure(1)
    pcolor(vec_temps, vec_freq, MAT_t_f_STFT_dB'); shading flat; caxis([30 90]); colorbar
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1?Pa2/Hz','interpreter','tex')
    set(gca,'FontSize',16)
    colormap jet
    %caxis([Lmin Lmax])
    
    if printFig ==true
        print([folderOut 'spectrogramChanelOne_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f1')
    end
end


% -------------------------------------------------------------------

if any( showFig == 2 )
    % azimut_NARW = atan2(y_NARW-y0,x_NARW-x0);
    % azimut_SHIP = atan2(y_SHIP-y0,x_SHIP-x0);
    
    figure(2)
    plot(vec_azimut*180/pi, Energie_dB,'k','LineWidth',2);
    hold on
    xlabel(' azimut (?)','interpret','tex')
    ylabel(' SPL dB re. 1?Pa','interpreter','tex')
    grid on
    set(gca,'FontSize',16)
    
    if printFig ==true
        print([folderOut 'energyFctAzimutXY_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f2')
    end
end

% -----------------------------------------------------------------
% Figure 3 : Energy fct azimut circle
if any( showFig == 3 )
    figure(3)
    R_source = 10;
    vec_azim = linspace(0,2*pi,1000);
    if strcmp(arr.arrOri , 'clock')
        x_cercle = R_source*sin(vec_azim);
        y_cercle = R_source*cos(vec_azim);
    else
        x_cercle = R_source*cos(vec_azim);
        y_cercle = R_source*sin(vec_azim);
    end
    
    plot(mean(xc),mean(yc),'k+')
    hold on
    plot(xc,yc,'r.','MarkerSize',12)
    plot(x_cercle,y_cercle,'k')
    for u = 1 : length(xc)
        text(xc(u)+0.3,yc(u),['h' num2str(u)])
    end
    
    xlabel(' x (m)')
    ylabel(' y (m)')
    grid on
    axis equal
    
    xxx = sin(vec_azimut).*Energie_NORM*R_source;
    yyy = cos(vec_azimut).*Energie_NORM*R_source;
    plot(xxx,yyy,'k','LineWidth',2)
    
    % Draw a north arrow
    quiver(9,7,0,3,'Color',[0 0 0],'LineWidth',1)
    text(9.2,9.5,'N')
    
    
    if printFig ==true
        print([folderOut 'energyFctAzimutCircle_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f3')
    end
    
    %{
    for u = 1 : length(xc)
        if strcmp(arr.arrOri, 'clock')
            a = atan2(xc(u),yc(u));
            xx = R_source*sin(a);
            yy = R_source*cos(a);
        else
            a = atan2(yc(u),xc(u));
            xx = R_source*cos(a);
            yy = R_source*sin(a);
        end
        plot(xx,yy,'r.');
        text(xx,yy,['h' num2str(u)])
    end
    %}
    % for u = 1 : length(x_NARW)
    %    text(x_NARW(u),y_NARW(u),['NARW' num2str(u)])
    % end
    % for u = 1 : length(x_SHIP)
    %    text(x_SHIP(u),y_SHIP(u),['SHIP' num2str(u)])
    % end
    
    
end % end showfig3

% ---------------------------------------------------------------------
% Figure 4: spectrogramme 1 voie
if any ( showFig == 4 )
    
    
    MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(indAziCible,:,:));
    for v = 1 : length(vec_f_INT)
        vec_pond = (MAT_POND_vs_h_freq(:,v)');
        vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
        vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
    end
    FFT_S_FV_int = zeros(1,floor(LFFT_FV/2)+1);
    FFT_S_FV_int(indi_f) = vec_sfft_FV;
    FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
    s_FV = ifft(FFT_S_FV);
    
    
    figure(4)
    t0=1;
    t =t0-Ns/2*1/fe+(0:Ns-1)*1/fe;
    subplot(5,1,1);
    plot(t,s_FV,'k')
    xlabel(' t (s)');
    %[vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(MAT_s_vs_t_h(:,1),0, fe, LFFT_spectro, REC, w_pond, fact_zp);
    [vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
    subplot(5,1,[2:5])
    pcolor(vec_temps_FV, vec_freq_FV,MAT_t_f_STFT_dB_FV'); shading flat; caxis([30 90]); colorbar
    ylim([fmin_int fmax_int])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(['VOIE n?' num2str(indAziCible) ', dB rel'])
    %caxis([Lmin Lmax])
    colormap jet
    
    if printFig ==true
        print([folderOut 'spectrgramoBestAngle_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f4')
    end
end





%  -------------------------------------------------------------------
% Figure 5 :  Spectogram N vois that cover 360 1 figure
if any ( showFig == 5)
    
    % Figure parameter
    nbSub = 4;
    sp.height=20; sp.width=40; sp.nbx=nbSub; sp.nby=nbSub;
    sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
    sp.spacex=0.5; sp.spacey=1;
    %sp.fracy = [0.1 0.3 0.3 0.3];
    sp.pos=subplot2(sp);
    
    % Needed some variable and matrix for the plot that are associated with
    % fig5
    
    
    % Reconstruct the signla for many direction
    Nvoie = nbSub^2;
    pas = ceil(length(vec_azimut)/nbSub^2);
    indNvoie = 1:pas:length(vec_azimut);
    
    % initialize
    vec_azim_nvoie = [];
    MAT_s_FV_vs_t_voie=zeros(Ns,length(1 : pas : length(vec_azimut)));
    
    
    for u = 1 : length(indNvoie)
        vec_azim_nvoie(u) = vec_azimut(indNvoie(u));
        MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(indNvoie(u),:,:));
        for v = 1 : length(vec_f_INT)
            vec_pond = (MAT_POND_vs_h_freq(:,v)');
            vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
            vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
        end
        FFT_S_FV_int = zeros(1,floor(LFFT_FV/2)+1);
        FFT_S_FV_int(indi_f) = vec_sfft_FV;
        FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
        s_FV = ifft(FFT_S_FV);
        MAT_s_FV_vs_t_voie(:,u)=s_FV;
    end
    
    % Call COMP_STFT once to have matrix size an vrialbe name
    t0=1;
    [vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
    
    
    % calcul des spectrogrammes des voies
    vec_freq_FV1 =(0:1:LFFT_spectro*fact_zp-1)*fe/(LFFT_spectro*fact_zp);
    indi_fV = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
    MAT_S_FV_vs_voie_t_f = zeros(Nvoie,length(vec_temps_FV),length(indi_fV));
    for u = 1 : Nvoie
        %     u
        %     Nvoie
        [vec_temps_FV, vec_freq_FV1, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(MAT_s_FV_vs_t_voie(:,u),t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
        indi_fV = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int))); %% THAT WERD MUST CHECK THIS OUT!!!!!!!
        MAT_S_FV_vs_voie_t_f(u,:,:)=MAT_t_f_STFT_dB_FV(:,indi_fV);
        %
    end
    vec_freq_FV = vec_freq_FV1(indi_fV);
    
    
    % Home bf
    
% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 70];           % [dB] C limite pcolor
spgm.im.dur = Nsample/10000;%'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 50;                   % [%] image window overlap
spgm.im.figvision = false  ;           % [true false] visiblity of figure before saveas jpf file
spgm.im.movm = 100;                  % Movmean parameter for the first panel
spgm.im.FontS = 14;                  % Image font Size
% spectrogram window parameters
spgm.win.dur = 0.18;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 4;                  % [s] spectrogram zero padding
spgm.fs = audioInfo.Fs;
spgm.im.ns = audioInfo.Ns;
%spgm.im.ns = fix(spgm.im.dur*acinfo.SampleRate);
%spgm.im.n = round(acinfo.TotalSamples/spgm.im.ns);
spgm.win.ns = fix(spgm.win.dur*spgm.fs);
spgm.win.wpond = kaiser(spgm.win.ns ,0.1102*(180-8.7));
spgm.win.val = spgm.win.wpond*sqrt(spgm.win.ns/sum(spgm.win.wpond.^2)); %w_pond    
spgm.win.novlp = fix(spgm.win.ovlp/100*spgm.win.ns);
spgm.win.nfft = max(256,2^nextpow2(spgm.win.ns*10));
spgm.im.fmin = 100;spgm.im.fmax = 200;

    [reconFFT, timeV, freqV] = beamForming('MLB', MAT_s_vs_t_h_INT , vec_azim_nvoie * 180 / pi, spgm);
    
    %compSpectro(MAT_S_FV_vs_voie_t_f(7,:,:),vec_freq_FV,vec_temps_FV,reconFFT(7,:,:),timeV,freqV,spgm)
    
    % ------ Figure plot --------
    cmapBW = flipud(gray(128));
    figure(5)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [sp.width sp.height]);
    set(gcf, 'PaperPositionMode', 'manual');
    
    meshNvoie = reshape(1:Nvoie,sp.nbx,sp.nby);
    
    for u = 1 :  Nvoie
        [sx,sy] = find(meshNvoie==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat;
        ylim([fmin_int fmax_int])
        
        if sx==1
            ylabel(' f (Hz)')
        else
            ax(u).YTickLabel = [];
        end
        if sy==sp.nby
            xlabel(' t (s)')
        else
            ax(u).XTickLabel = []
        end
        
        title(['V ' num2str(u) 'Az (?) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
        caxis([Lmin Lmax])
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{iFile})
    
    
    if printFig ==true
        print([folderOut 'subplotSpectroNVoie_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f5')
    end
end



%  -------------------------------------------------------------------
% Figure 6 :  Spectogram N vois over cible 1 figure
if any ( showFig == 6)
    % Figure parameter
    nbSub = 4;
    Nvoie = nbSub^2 ;
    sp.height=20; sp.width=40; sp.nbx=nbSub; sp.nby=nbSub;
    sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
    sp.spacex=0.5; sp.spacey=1;
    %sp.fracy = [0.1 0.3 0.3 0.3];
    sp.pos=subplot2(sp);
    
    % Needed some variable and matrix for the plot that are associated with
    % fig5
    % Reconstruct the signla for many direction
    pas = 2;
    indNVoie = indAziCible-(Nvoie/2)*pas+1:pas:indAziCible+(Nvoie/2)*pas;
    
    % initialize
    vec_azim_nvoie = [];
    MAT_s_FV_vs_t_voie=zeros(Ns,length(1 : pas : length(vec_azimut)));
    
    
    for u = 1 : length(indNVoie)
        vec_azim_nvoie(u) = vec_azimut(indNVoie(u));
        MAT_POND_vs_h_freq =squeeze(MAT_POND_vs_azim_h_freq(indNVoie(u),:,:));
        for v = 1 : length(vec_f_INT)
            vec_pond = (MAT_POND_vs_h_freq(:,v)');
            vec_sfft  = MAT_ffts_vs_f_h_INT_INT(v,:);
            vec_sfft_FV(v) = sum(vec_pond.*vec_sfft);
        end
        FFT_S_FV_int = zeros(1,floor(LFFT_FV/2)+1);
        FFT_S_FV_int(indi_f) = vec_sfft_FV;
        FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
        s_FV = ifft(FFT_S_FV);
        MAT_s_FV_vs_t_voie(:,u)=s_FV;
    end
    
    % Call COMP_STFT once to have matrix size an vrialbe name
    t0=1;
    [vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
    
    
    % calcul des spectrogrammes des voies
    vec_freq_FV1 =(0:1:LFFT_spectro*fact_zp-1)*fe/(LFFT_spectro*fact_zp);
    indi_fV = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
    MAT_S_FV_vs_voie_t_f = zeros(Nvoie,length(vec_temps_FV),length(indi_fV));
    for u = 1 : Nvoie
        %     u
        %     Nvoie
        [vec_temps_FV, vec_freq_FV1, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(MAT_s_FV_vs_t_voie(:,u),t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);
        indi_fV = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
        MAT_S_FV_vs_voie_t_f(u,:,:)=MAT_t_f_STFT_dB_FV(:,indi_fV);
        %
    end
    vec_freq_FV = vec_freq_FV1(indi_fV);
    
    
    
    
    % ------ Figure plot --------
    cmapBW = flipud(gray(128));
    figure(6)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [sp.width sp.height]);
    set(gcf, 'PaperPositionMode', 'manual');
    
    meshNvoie = reshape(1:Nvoie,sp.nbx,sp.nby);
    
    for u = 1 :  Nvoie
        [sx,sy] = find(meshNvoie==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat;
        ylim([fmin_int fmax_int])
        
        if sx==1
            ylabel(' f (Hz)')
        else
            ax(u).YTickLabel = [];
        end
        if sy==sp.nby
            xlabel(' t (s)')
        else
            ax(u).XTickLabel = []
        end
        
        title(['V ' num2str(u) 'Az (?) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
        caxis([Lmin Lmax])
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{iFile})
    
    
    
    if printFig ==true
        print([folderOut 'subplotSpectroNVoieCible_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f6')
    end
end




% ---------------------------------------------------------------
% Figure 9 : spectogram de N voies to video
if any (showFig == 99 )
    % This need to be recode !!!! ----------------------
    
    
    spie = 1;
    pas = floor(length(vec_azimut)/NB_voies_visees);
    MAT_s_FV_vs_t_voie=zeros(Ns,length(1 : pas : length(vec_azimut)));
    
    for u = 1 : pas : length(vec_azimut)
        
        vec_azim_nvoie(spie) = vec_azimut(u);
        indice_voie(spie) = u;
        
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
        
        MAT_s_FV_vs_t_voie(:,spie)=s_FV;
        spie = spie + 1;
    end
    
    
    % calcul des spectrogrammes des voies
    Nvoie = length(vec_azim_nvoie);
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
    
    
    % trac? des spectrogrammes
    
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
        %      Nvoie
        figure(9)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat; caxis([30 90])
        ylim([fmin_int fmax_int])
        xlabel(' t (s)')
        ylabel(' f (Hz)')
        title(['V ' num2str(u) 'Az (?) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
        caxis([Lmin Lmax])
        colorbar
        colormap jet
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
        
        if printFig ==true
            if ~isfolder([folderOut 'videoSnap\']); disp(['Creating output folder: ' folderOut 'videoSnap\']); mkdir([folderOut 'videoSnap\']); end
            print([folderOut 'videoSnap\spectrogramNAngle_' outName '_' num2str(vec_azim_nvoie(u)* 360 / 2/pi)   'd.png'], '-r150','-dpng')
        end
        
        
        h = gcf;
        F = getframe(h);
        writeVideo(writerObj,F)
        close
        
    end
    
    close(writerObj)
    
end % end showfig 5
% -----------------------------------------------------------------------