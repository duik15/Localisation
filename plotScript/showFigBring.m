% This is just a script containing the figure for each timestepd on locateBring 
% ---------------------------------------------------------------------
% Figure 1 : spectogram
% Figure 2 : Energy fct azimut
% Figure 3 : Energy fct azimut circle
% Figure 4: spectrogramme 1 voie
% Figure 5 :  Spectogram N vois 1 figure
% Figure 9 : spectogram de N voies in a forloop for video


% ------------------------ Figure in loop -----------------------------

% Figure 1 : spectogram
if any( showFig == 1 )
    figure(1)
    pcolor(vec_temps, vec_freq, MAT_t_f_STFT_dB'); shading flat; caxis([30 90]); colorbar
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1µPa2/Hz','interpreter','tex')
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
    xlabel(' azimut (°)','interpret','tex')
    ylabel(' SPL dB re. 1µPa','interpreter','tex')
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
    R_source = 5000;
    vec_azim = linspace(0,2*pi,1000);
    if strcmp(arrOri , 'clock')
        x_cercle = R_source*sin(vec_azim);
        y_cercle = R_source*cos(vec_azim);
    else
        x_cercle = R_source*cos(vec_azim);
        y_cercle = R_source*sin(vec_azim);
    end
    
    plot(x0,y0,'k+')
    hold on
    plot(xc,yc,'r.','MarkerSize',12)
    plot(x_cercle,y_cercle,'k')
    
    
    xlabel(' x (m)')
    ylabel(' y (m)')
    grid on
    axis equal
    for u = 1 : length(xc)
        if strcmp(arrOri, 'clock')
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
    % for u = 1 : length(x_NARW)
    %    text(x_NARW(u),y_NARW(u),['NARW' num2str(u)])
    % end
    % for u = 1 : length(x_SHIP)
    %    text(x_SHIP(u),y_SHIP(u),['SHIP' num2str(u)])
    % end
    
    
    xxx = sin(vec_azimut).*Energie_NORM*R_source;
    yyy = cos(vec_azimut).*Energie_NORM*R_source;
    plot(xxx,yyy,'k','LineWidth',2)
    
    if printFig ==true
        print([folderOut 'energyFctAzimutCircle_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f3')
    end
end % end showfig3

% ---------------------------------------------------------------------
% Figure 4: spectrogramme 1 voie
if any ( showFig == 4 )
    
    % formation d'une voie
    [~, iamax ]  =max(Energie_dB);
    azimut_cible = vec_azimut(iamax);
    [ val indi_azimut_vise] = min(abs(vec_azimut - azimut_cible));

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
    
    if printFig ==true
        print([folderOut 'spectrgramoBestAngle_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f4')
    end
end





%  -------------------------------------------------------------------
% Figure 5 :  Spectogram N vois 1 figure
if any ( showFig == 5)
    
    % Needed some variable and matrix for the plot that are associated with fig9
    if ~exist('MAT_s_FV_vs_t_voie') 
        
        
       

        
        spie = 1;
        nbSub = 4;
        pas = ceil(length(vec_azimut)/nbSub^2);
        MAT_s_FV_vs_t_voie=zeros(Ns,length(1 : pas : length(vec_azimut)));
        
        for u = 1 : pas : length(vec_azimut)
            
            vec_azim_voie(spie) = vec_azimut(u);
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
       
       % Call COMP_STFT once to have matrix size an vrialbe name 
       t0=1;
       [vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, LFFT_spectro, REC, w_pond, fact_zp);

        
        % calcul des spectrogrammes des voies
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
        
        
    end
    
    % ------ Figure plot ---------
    figure(5)
    for u = 1 : 1: Nvoie
        %      u
        %      Nvoie
        subplot(nbSub,nbSub,u)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat;
        ylim([fmin_int fmax_int])
        xlabel(' t (s)')
        ylabel(' f (Hz)')
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(180/pi*vec_azim_voie(u)))])
        caxis([Lmin Lmax])
        colormap jet
    end
    
    if printFig ==true
        print([folderOut 'subplotSpectroNAngle_' outName '_p' num2str(iFile)  '.png'], '-r150','-dpng', '-f5')
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
        
        vec_azim_voie(spie) = vec_azimut(u);
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
   
    
    % tracé des spectrogrammes

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
        
        if printFig ==true
        if ~isfolder([folderOut 'videoSnap\']); disp(['Creating output folder: ' folderOut 'videoSnap\']); mkdir([folderOut 'videoSnap\']); end
        print([folderOut 'videoSnap\spectrogramNAngle_' outName '_' num2str(vec_azim_voie(u)* 360 / 2/pi)   'd.png'], '-r150','-dpng', '-f9')
        end
        
        
        h = gcf;
        F = getframe(h);
        writeVideo(writerObj,F)
        close
        
    end
    
    close(writerObj)
   
end % end showfig 5
% -----------------------------------------------------------------------