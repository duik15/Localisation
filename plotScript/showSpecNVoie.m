    
    % Number of direction
    nbV =12;
    indV = ceil(length(vec_azimut)/nbV) * (1:nbV);
    azimutV = vec_azimut(indV) * 180 / pi; 
    % Figure parameter
    sp.height=20; sp.width=40; sp.nbx=2; sp.nby=nbV/sp.nbx;
    sp.ledge=5; sp.redge=2; sp.tedge=2; sp.bedge=2;
    sp.spacex=0.3; sp.spacey=0.2;
    %sp.fracy = [0.1 0.3 0.3 0.3];
    sp.pos=subplot2(sp);
    %%
    [reconFFT, timeV, freqV] = reconstruct(nbV,MAT_POND_vs_azim_h_freq, MAT_ffts_vs_f_h_INT_INT, vec_azimut, vec_f_INT, spec);
    
    
    %%
    clf
    % ------ Figure plot --------
    figure(1)
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [sp.width sp.height]);
    set(gcf, 'PaperPositionMode', 'manual');
    
  meshNvoie = reshape(1:nbV,sp.nby,sp.nbx)


    
    for u = 1 :  nbV
        [sy,sx] = find(meshNvoie==u);
        ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
        %subplot(nbSub,nbSub,u)
        %M =  (squeeze(recon(u,:,:)))';
        pcolor(timeV, freqV,(squeeze(reconFFT(u,:,:)))'); shading flat;
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
        
        if u <= nbV/2
           tt = title([ num2str(azimutV(u)+1) '$^{\circ}$'],'interpreter','latex');
           tt.Position = [ timeV(1)-0.2*(timeV(end)-timeV(1))  (round(freqV(end)-freqV(1)))/2+freqV(1) ]
        else
           tt = title([ num2str(azimutV(u)+1) '$^{\circ}$'],'interpreter','latex');
           tt.Position = [ timeV(end)+0.07*(timeV(end)-timeV(1))  (round(freqV(end)-freqV(1)))/2+freqV(1)-5 ]
        end
        
        caxis([spec.Lmin spec.Lmax])
        colormap jet
        %colormap(cmapBW)
        grid on
        set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
    end
    
    sgtitle(fileList{iFile})
    
    
    if printFig ==true
        print([folderOut 'subplotFormationVoie_360_' outName '_p' num2str(iFile)  '.png'], '-r300','-dpng', '-f1')
    end