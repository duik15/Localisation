
% pcolormap of the matrice energie as function of time and azimute
% Data p,g and result are plot also

figure(2)
ax=gca;
colorO = ax.ColorOrder;
legGroup = [];
legName = {};

pcolor(arrIncol(ptime), 1:360, 10*log10(matEnergie'))
shading flat
cb = colorbar('location','eastoutside');
ylabel(cb, 'Energy (dB)');
colormap jet

% Ping points
if exist('p')
hold on
hpm = plot(p.ptime,p.angleM,'p','color','r','markersize',15,'markerfacecolor','r');
legGroup = [legGroup hpm];
legName = {legName{:}, 'Ping azimuth'};
end

% GPS tracks points
if exist('g')
hold on
hgr = plot(g.time,g.angleR,'o','color','k','markersize',5,'markerfacecolor','k');
legGroup = [legGroup hgr];
legName = {legName{:}, 'Real boat azimuht'};
end

% Result points
hm = plot(ptime,angleM,'o','color','r','markersize',5,'markerfacecolor','r');
legGroup = [legGroup hm];
legName = {legName{:} 'Boat azimuth'};


% Axes paremeter
ylabel('Azimuth $^{\circ}$','interpreter','latex')
xlabel('Time')
    
    
leg = legend(legGroup,legName);
%leg.Box ='off';

print('-dpng','-r150',[folderOut 'matEnergieNAzigrame_' outName '.png' ])
