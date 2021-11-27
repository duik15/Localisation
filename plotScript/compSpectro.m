%function compSpectro(Pdb1,t1,f1,Pdb2,t2,f2,spgm)
% Debug 
Pdb1 = squeeze(MAT_S_FV_vs_voie_t_f(7,:,:))';
t1 =vec_temps_FV;
f1=vec_freq_FV; 
Pdb2 = fliplr(squeeze(reconFFT(7,:,:))');
t2 = timeV;
f2= freqV;

close all


%Pdb1 = squeeze(Pdb1);
%Pdb2 = squeeze(Pdb2)';

% Compare two spectrogram

sp2.height=20; sp2.width=40; sp2.nbx=1; sp2.nby=2;
sp2.ledge=3; sp2.redge=3.5; sp2.tedge=2; sp2.bedge=2;
sp2.sp2acex=0.3; sp2.sp2acey=0.3;
%sp2.fracy = [0.1 0.3 0.3 0.3];
sp2.pos=subplot2(sp2);

figure(1)
sizeFont= 24;
set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen

% Creat meash grid
meshNvoie = reshape(1:2,2,1);

for u=1:2
    [sy,sx] = find(meshNvoie==u);
    ax(u)=axes('position',sp2.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
    %ax = subplot(2,1,1)
    %subplot(nbSub,nbSub,u)
    %M =  (squeeze(recon(u,:,:)))';
    
    %Plot
    if u==1
    pcolor(t1, f1,Pdb1); shading flat;
    else
        pcolor(t2, f2,Pdb2); shading flat;
    end
    spgm.im.fmin
    spgm.im.fmax
    ylim([spgm.im.fmin spgm.im.fmax])
    set(ax(u),'FontSize',14)
    if sx==1
        ylabel(' f (Hz)','FontSize',sizeFont)
        %set(gca,'Xtick',[150 160 170 180 190],'xtickLabel',num2cell([150 160 170 180 190]))
    else
        ax(u).YTickLabel = [];
    end
    if sy==2
        xlabel(' t (s)','FontSize',sizeFont)
    else
        ax(u).XTickLabel = [];
    end
    
    %if u <= nbV/2
        %tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        %tt.Position = [ timeWin(1)-0.12*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1) ];
    %else
        %tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        %tt.Position = [ timeWin(end)+0.15*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1)-5 ];
    %end
    
    caxis(spgm.im.clims)
    colormap jet
    %colormap(cmapBW)
    grid on
    
end  

%sgtitle([file{i} '  |  ' num2str(floor(timeWin(1))) ' s to ' num2str(ceil(timeWin(end))) ' s  |  File ' num2str(j) '/' num2str(spec.nbIm) ] ,'FontSize', sizeFont+4)

% Color bar
cbPos = [sp2.pos{end}(1)+sp2.pos{end}(3)+0.01 sp2.pos{end}(2) 0.01 (sp2.height - sp2.bedge - sp2.tedge)/sp2.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
ax(end).Position = sp2.pos{end};
ylabel(cb,'Power (dB)')

%print('-dpng','-r150',[folderOut 'debug_BF_' file{i}(1:end-4) '_s' num2str(floor(timeWin(1))) '_e' num2str(ceil(timeWin(end))) '.png' ])