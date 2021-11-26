

% ------ Figure plot --------
figure(1)
sizeFont= 24;
set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen

% Creat meash grid
meshNvoie = reshape(1:nbV,sp.nby,sp.nbx);

for u = 1 :  nbV
    [sy,sx] = find(meshNvoie==u);
    ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
    %subplot(nbSub,nbSub,u)
    %M =  (squeeze(recon(u,:,:)))';
    pcolor(timeV + timeWin(1), freqV,(squeeze(reconFFT(u,:,:)))'); shading flat;
    ylim([spec.fmin spec.fmax])
    set(ax(u),'FontSize',14)
    if sx==1
        ylabel(' f (Hz)','FontSize',sizeFont)
        %set(gca,'Xtick',[150 160 170 180 190],'xtickLabel',num2cell([150 160 170 180 190]))
    else
        ax(u).YTickLabel = [];
    end
    if sy==sp.nby
        xlabel(' t (s)','FontSize',sizeFont)
    else
        ax(u).XTickLabel = [];
    end
    
    if u <= nbV/2
        tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        tt.Position = [ timeWin(1)-0.12*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1) ];
    else
        tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        tt.Position = [ timeWin(end)+0.15*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1)-5 ];
    end
    
    caxis([spec.Lmin spec.Lmax])
    colormap jet
    %colormap(cmapBW)
    grid on
    set(gca,'layer','top','gridAlpha',0.5,'XMinorGrid','on','YMinorGrid','on','MinorGridAlpha',0.3,'XMinorTick','on','YMinorTick','on')
end

% Title
sgtitle([file{i} '  |  ' num2str(floor(timeWin(1))) ' s to ' num2str(ceil(timeWin(end))) ' s  |  File ' num2str(j) '/' num2str(spec.nbIm) ] ,'FontSize', sizeFont+4)


% Color bar
cbPos = [sp.pos{end}(1)+sp.pos{end}(3)+0.01 sp.pos{end}(2) 0.01 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
ax(end).Position = sp.pos{end};
ylabel(cb,'Power (dB)')

print('-dpng','-r150',[folderOut 'BF_' file{i}(1:end-4) '_s' num2str(floor(timeWin(1))) '_e' num2str(ceil(timeWin(end))) '.png' ])

