    % pcolormap of the matrice energie
figure(2)

ax=gca;
pcolor(ptime, 1:360, 10*log10(matEnergie'))
shading flat
cb = colorbar('location','eastoutside');
ylabel(cb, 'Energy (dB)');
colormap jet

hold on
hm = plot(ptime,angleM,'o','color','k','markersize',5);

% Axes paremeter
ylabel('Azimut $^{\circ}$','interpreter','latex')
xlabel('Time')

leg = legend([hm],'First  Lobe','Location','northwest');
%leg.Box ='off';

print('-dpng','-r150',[folderOut 'pcolorAngleTime_' outName '.png' ])