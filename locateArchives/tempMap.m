% This script plot arrival angle modify with sidelobe


clf
% Linéaire regression
x = 1:10;
p = polyfit(x,angleM,1);
angleFit = polyval(p,x);  


figure(7)
hA = plot(x,angleM,'-.ok');
hold on
hTh = plot(x,0:36:360-36,'k-');
hFit = plot(x,angleFit);

% Plot some side lobe
for ii=2:3
   hSl(ii-1) = plot(x,angleA(:,ii),'o') 
end

xlabel('Ping number')
ylabel('Angle of arrival')
ax = gca;
ax.XTickLabel = pNb;
leg= legend([hA,hTh,hFit, hSl(1), hSl(2)],{'Ping','Theoritical curve','Lin. Reg. Curve','Second lobe', 'Third Lobe'},'Location','northwest')

print([folderOut 'angleOfArrival_2N3Lobe_' outName '.png'],'-r150','-dpng')


%%

% This script plot arrival angle modify with sidelobe
lobeNb = [1 1 1 3 1 1 2 3 1 2];
for ii= 1:nbF
    angleBest(ii) = angleA(ii,lobeNb(ii)); 
end
clf
% Linéaire regression
x = 1:10;
p = polyfit(x,angleM,1);
angleFit = polyval(p,x);  


figure(7)
hA = plot(x,angleM,'-.ok');
hold on
hTh = plot(x,0:36:360-36,'k-');
hFit = plot(x,angleFit);

% Plot some side lobe
for ii=2:3
   hSl(ii-1) = plot(x,angleA(:,ii),'o') 
end

% Best angle
hBest = plot(x,angleBest,'x');
p = polyfit(x,angleBest,1);
angleFit = polyval(p,x);  
hBFit = plot(x,angleFit);

xlabel('Ping number')
ylabel('Angle of arrival')
ax = gca;
ax.XTickLabel = pNb;
leg= legend([hA,hTh,hFit, hSl(1), hSl(2),hBest,hBFit],{'Ping','Theoritical curve','Lin. Reg. Curve','Second lobe', 'Third Lobe','Best lobe','Best lobe Fit'},'Location','northwest')

print([folderOut 'angleOfArrival_bestLobe_' outName '.png'],'-r150','-dpng')
