%% ---------------- Global figures -------------------------------


if any ( showFig == 7)

if sum(diff(angleM(3:end-2)))< 0
    pNb = 10:-1:1;
    angleA = fliplr(aMax3);
else
    pNb = 1:10;
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

end

if any ( showFig ==8)
figure(8)

plot(xc,yc,'o')


for ic = 1:length(xc)
text(xc(ic)+0.2,yc(ic),num2str(ic),'fontsize',14)
end
end

if any ( showFig == 98)
    
    
    figure(9)
    
    plot(angleM,'--ok')
    xlabel('Ping number')
    ylabel('Angle of arrival')
    print([folderOut 'angleOfArrivalbad_' outName '.png'],'-r150','-dpng')
    
    
    
     [asort, isort] = sort(angleM);
    
    figure(7)
    
    plot(asort,'--ok')
    hold on
    plot([])
    xlabel('Ping number')
    ylabel('Angle of arrival')
    
    % Thick sort
    ax=gca;
    
    vec10 = 1:10;
    pingSort = vec10(isort)
    pingSortStr = strsplit(num2str(vec10(isort)));
    ax.XTickLabel = pingSortStr;
    print([folderOut 'angleOfArrival_sort_' outName '.png'],'-r150','-dpng')
end