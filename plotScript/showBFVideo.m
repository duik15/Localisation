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
        FFT_S_FV_int = zeros(1,floor(spec.Ns/2)+1);
        FFT_S_FV_int(indi_f) = vec_sfft_FV;
        FFT_S_FV = [ FFT_S_FV_int conj(fliplr(FFT_S_FV_int(2:length(FFT_S_FV_int)-1)))];
        s_FV = ifft(FFT_S_FV);
        
        MAT_s_FV_vs_t_voie(:,spie)=s_FV;
        spie = spie + 1;
    end
    
    
    % calcul des spectrogrammes des voies
    Nvoie = length(vec_azim_nvoie);
    % ATTENITON À CETTE LIGNES!!!
    vec_freq_FV1 =(0:1:spec.winSz*spec.zp-1)*fe/(spec.winSz*spec.zp);
    indi_f = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
    MAT_S_FV_vs_voie_t_f = zeros(Nvoie,length(vec_temps_FV),length(indi_f));
    for u = 1 : Nvoie
        %     u
        %     Nvoie
        [vec_temps_FV, vec_freq_FV1, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(MAT_s_FV_vs_t_voie(:,u),t0-Ns/2*1/fe, fe, spec.winSz, spec.rec, spec.wpond, spec.zp);
        indi_f = find( ((vec_freq_FV1>= fmin_int)&(vec_freq_FV1 <= fmax_int)));
        if u==1
           saveindi = indi_f; 
        end
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
        figure(9)
        M =  squeeze(MAT_S_FV_vs_voie_t_f(u,:,:));
        pcolor(vec_temps_FV, vec_freq_FV,M'); shading flat; caxis([30 90])
        ylim([fmin_int fmax_int])
        xlabel(' t (s)')
        ylabel(' f (Hz)')
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(180/pi*vec_azim_nvoie(u)))])
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