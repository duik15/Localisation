
if sum(diff(angleM(3:end-2)))< 0
    pNb = nbF:-1:1;
else
    pNb = 1:nbF;
end    
    
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
for ii=2:2
   hSl = plot(x,angleA(:,ii),'o') 
end

xlabel('Ping number')
ylabel('Angle of arrival')
ax = gca;
ax.XTickLabel = pNb;
leg= legend([hA,hTh,hFit, hSl],{'Ping','Theoritical curve','Lin. Reg. Curve','Secondary lobe'},'Location','northwest')

print([folderOut 'angleOfArrival_' outName '.png'],'-r150','-dpng')