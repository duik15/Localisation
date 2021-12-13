% Show real quick the recording
   
i_ch = 1; 

sp.height=30; sp.width=30; sp.nbx=1; sp.nby=2;
    sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
    sp.spacex=0.5; sp.spacey=1;
    sp.pos=subplot2(sp);

% Get the spectogram of first chanel
    [~,freq1,time1,tmp] = spectrogram(wav.pa(:,i_ch),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
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
        print([folderOut 'spectrogramChanel_' num2str(i_ch)  '_' outName '_p' num2str(ifile)  '.png'], '-r150','-dpng', '-f1')
    end
    
