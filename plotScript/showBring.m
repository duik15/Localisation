% This is just a script containing the figure for each timestepd on locateBring
% ---------------------------------------------------------------------
% Figure 1 : spectogram
% Figure 2 : Energy fct azimut
% Figure 3 : Energy fct azimut circle
% Figure 4: spectrogramme 1 voie
% Figure 5 :  Spectogram N vois that cover 360deg 1 figure
% Figure 6 :  Spectogram N vois around azi cible 1 figure

% ------------------- General Parameter ------------------------------
% formation de voie
if ~exist('aziCible') || strcmp(aziCible,'max')
    aziCible = angleM;
end
[ ~, indAziCible] = min(abs(azimut360 - aziCible));
aziCible2 =  aziCible();
arrDiameter = 10;

% ------------------------ Figure in loop -----------------------------
% -----------------------------------------------------------------
% Figure 1 : spectogram
% -----------------------------------------------------------------
if any( showFig == 1 )
    
    % Get the spectogram of first chanel
    [~,freq1,time1,tmp] = spectrogram(wav.db(:,1),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    Pdb1 = 10*log10(tmp);
    
    figure(1)
    pcolor(time1, freq1, Pdb1); shading flat;
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1µPa2/Hz','interpreter','tex')
    set(gca,'FontSize',16)
    colormap jet
    caxis(spgm.im.clims)
    cb = colorbar;
    ylabel(cb,'Power (dB)')
    if printFig ==true
        print([folderOut 'spectrogramChanelOne_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f1')
    end
end



% -------------------------------------------------------------------
% Figure 2 : Energy fct azimut
% -------------------------------------------------------------------
if any( showFig == 2 )
    
    figure(2)
    plot(azimut360, Energie_dB,'k','LineWidth',2);
    hold on
    xlabel(' azimut (°)','interpret','tex')
    ylabel(' SPL dB re. 1µPa','interpreter','tex')
    grid on
    set(gca,'FontSize',16)
    xlim([0 360])
    if printFig ==true
        print([folderOut 'energyFctAzimutXY_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f2')
    end
end



% -----------------------------------------------------------------
% Figure 3 : Energy fct azimut circle
% -----------------------------------------------------------------
if any( showFig == 3 )
    
    aziCircle = linspace(0,2*pi,1000);
    
    figure(3)
    plot(mean(arr.xh),mean(arr.yh),'k+')
    hold on
    plot(arr.xh,arr.yh,'r.','MarkerSize',12)
    plot(arrDiameter*sin(aziCircle),arrDiameter*cos(aziCircle),'k')
    for u = 1 : length(arr.xh)
        text(arr.xh(u)+0.3,arr.yh(u),['h' num2str(u)])
    end
    
    xlabel(' x (m)')
    ylabel(' y (m)')
    grid on
    axis equal
    
    xxx = sin(azimut360* pi /180).*Energie_NORM*arrDiameter;
    yyy = cos(azimut360* pi /180).*Energie_NORM*arrDiameter;
    plot(xxx,yyy,'k','LineWidth',2)
    
    % Draw a north arrow
    quiver(9,7,0,3,'Color',[0 0 0],'LineWidth',1)
    text(9.2,9.5,'N')
    
    
    if printFig ==true
        print([folderOut 'energyFctAzimutCircle_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f3')
    end
    
end % end showfig3




% ---------------------------------------------------------------------
% Figure 4: spectrogramme 1 voie
% ---------------------------------------------------------------------
if any ( showFig == 4 )
    clear ax
    %[vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(MAT_s_vs_t_h(:,1),0, fe, spec.winSz, spec.rec, spec.wpond, spec.zp);
    %[vec_temps_FV, vec_freq_FV, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB_FV] = COMP_STFT_snapshot(s_FV,t0-Ns/2*1/fe, fe, spec.winSz, spec.rec, spec.wpond, spec.zp);
    
    % Get spectro of cible
    [PdbC, timeC, freqC,reconC] = beamForming(arrID, wav.db , angleM , spgm,'specmethod','spectro');
    PdbC = squeeze(PdbC);
    
    disp('Wong time in figure4')
    %tRecon =1-spgm.fs/2*1/spgm.fs+(0:spgm.ns-1)*1/spgm.fs;
    t1=floor(spgm.ns/4);
    tRecon =[time_s(t1:end) time_s(1:t1-1) ];
    figure(4)
    ax(1) = subplot(5,1,1);
    plot(tRecon,reconC,'k')
    title(['Azimut ' num2str(azimut360(indAziCible)) '°'])
    
    ax(2) = subplot(5,1,[2:5]);
    pcolor(timeC, freqC,PdbC'); shading flat;
    ylim([spgm.im.freqlims])
    xlabel(' Time (s)')
    ylabel(' Frequency (Hz)')
    caxis(spgm.im.clims)
    
    %caxis([Lmin Lmax])
    colormap jet
    
    % Coordinate axis
    ax(1).XLim = ax(2).XLim;
    
    % Colorbar
    tmppos = ax(2).Position;
    cb = colorbar;
    ax(2).Position = tmppos;
    ylabel(cb,'Power (dB)')
    
    
    
    if printFig ==true
        print([folderOut 'spectrgramoBestAngle_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f4')
    end
end





%  -------------------------------------------------------------------
% Figure 5 :  Spectogram N vois that cover 360 1 figure
% -----------------------------------------------------------------
if any ( showFig == 5)
    
    % Figure parameter
    nbSub = 4;
    Nvoie =  nbSub^2;
    sp.height=30; sp.width=30; sp.nbx=nbSub; sp.nby=nbSub;
    sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
    sp.spacex=0.5; sp.spacey=1;
    sp.pos=subplot2(sp);
    
    % Azimut to plot
    azimutV = 0:360/((nbSub^2)-1):360;
    
    % Get the beamforming spectrogram
    [Pdb, timeV, freqV] = beamForming(arrID, wav.db , azimutV, spgm,'specmethod','spectro');
    
    
    % ------ Figure plot --------
    figure(5)
    set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen
    meshNvoie = reshape(1:Nvoie,sp.nbx,sp.nby);
    
    for u = 1 :  Nvoie
        [sx,sy] = find(meshNvoie==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        M =  squeeze(Pdb(u,:,:));
        pcolor(timeV, freqV,M'); shading flat;
        ylim([spgm.im.freqlims])
        
        if sx==1
            ylabel(' f (Hz)')
        else
            ax(u).YTickLabel = [];
        end
        if sy==sp.nby
            xlabel(' t (s)')
        else
            ax(u).XTickLabel = [];
        end
        
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(azimutV(u)))])
        caxis(spgm.im.clims)
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{ifile},'interpreter','none')
    
    % Color bar
    cbPos = [sp.pos{end}(1)+sp.pos{end}(3)+0.01 sp.pos{end}(2) 0.01 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
    cb= colorbar;%('Position',cbPos);
    cb.Position = cbPos;
    ax(end).Position = sp.pos{end};
    ylabel(cb,'Power (dB)')
    
    if printFig ==true
        print([folderOut 'subplotSpectroNVoie_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f5')
    end
end



%  -------------------------------------------------------------------
% Figure 6 :  Spectogram N vois over cible 1 figure
if any ( showFig == 6)
    
    
    % Figure parameter
    nbSub = 4;
    Nvoie =  nbSub^2;
    sp.height=30; sp.width=30; sp.nbx=nbSub; sp.nby=nbSub;
    sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
    sp.spacex=0.5; sp.spacey=1;
    sp.pos=subplot2(sp);
    
    % Azimut para
    deltaA = 2;
    azimutV = angleM - deltaA * Nvoie/2:deltaA:angleM + deltaA * Nvoie/2;
    
    % Get the beamforming spectrogram
    [Pdb, timeV, freqV] = beamForming(arrID, wav.db , azimutV, spgm,'specmethod','spectro');
    
    %[recon, timeV, freqV] = beamForming(matPondahf, matFFT, azimut, freq, spec);
    
    % ------ Figure plot --------
    figure(6)
    set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen
    meshNvoie = reshape(1:Nvoie,sp.nbx,sp.nby);
    
    for u = 1 :  Nvoie
        [sx,sy] = find(meshNvoie==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        M =  squeeze(Pdb(u,:,:));
        pcolor(timeV, freqV,(M')); shading flat;
        ylim([spgm.im.freqlims])
        
        if sx==1
            ylabel(' f (Hz)')
        else
            ax(u).YTickLabel = [];
        end
        if sy==sp.nby
            xlabel(' t (s)')
        else
            ax(u).XTickLabel = [];
        end
        
        title(['V ' num2str(u) 'Az (°) ' num2str(floor(azimutV(u)))])
        caxis(spgm.im.clims)
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{ifile},'interpreter','none')
    
    % Color bar
    cbPos = [sp.pos{end}(1)+sp.pos{end}(3)+0.01 sp.pos{end}(2) 0.01 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
    cb= colorbar;%('Position',cbPos);
    cb.Position = cbPos;
    ax(end).Position = sp.pos{end};
    ylabel(cb,'Power (dB)')
    
    
    if printFig ==true
        print([folderOut 'subplotSpectroNVoieCible_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f6')
    end
end




% ---------------------------------------------------------------
