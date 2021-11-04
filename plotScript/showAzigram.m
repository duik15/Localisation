% This figure shwo the angle of detection throught time
% Azigrame/ 


clf

ploth=15;
plotw=30;

figure(1)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotw ploth]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotw ploth]);

% Get axis properties
ax=gca;
colorO = ax.ColorOrder;
legGroup = [];
legName = {};

if exist('g') % GPS track
    hr = plot(g.time,g.angleR,'o','color','k','markerfacecolor','k','markersize',5);
    legGroup = [legGroup hr];
    legName = {legName{:}, 'True azimuth'};
elseif exist('angleR')
    hr = plot(ptime,angleR,'o','color','k','markerfacecolor','k','markersize',5);
    legGroup = [legGroup hr];
    legName = {legName{:}, 'True azimuth'};
end
    
hold on
hm = plot(ptime,angleM,'o','color',colorO(2,:),'markerfacecolor',colorO(2,:),'markersize',3);
h2 = plot(ptime,angleA(:,2),'.','color',colorO(3,:),'markersize',6);
h3 = plot(ptime,angleA(:,3),'.','color',colorO(4,:),'markersize',6);
legGroup = [legGroup hm h2 h3];
legName = {legName{:}, 'First lobe','Second lobe','Third lobe'};

if exist('p') % Ping emision
hpm = plot(p.ptime,p.angleM,'p','color','r','markersize',20,'markerfacecolor','r');
legGroup = [legGroup hpm];
legName = {legName{:}, 'Located ping'};
end

ylabel('Azimuth $^\circ$','interpreter','latex')
xlabel('Time')

if exist(limBol) && limBol == true
    xlim(xlimit)
end

leg = legend(legGroup,legName,'Location','northwest');
print('-dpng','-r150',[folderOut 'angleTime_' outName '.png' ])

