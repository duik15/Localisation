function compSpectro(Pdb1,t1,f1,Pdb2,t2,f2,spgm, varargin)

while ~isempty(varargin)
        switch lower(varargin{1})
            case 'time'
                time = varargin{2};
            case 'type'
                type = varargin{2};
            case 'typeplot'
                 typePlot = varargin{2};
            case 'name'
                nameSub = varargin{2};
            case 'outname'
                outName = varargin{2};
            otherwise
                error(['Can''t understand property: ' varargin{1}])
        end
            varargin(1:2)=[];
end

% Debug
%Pdb1 = squeeze(Pdb1);
%Pdb2 = squeeze(Pdb2)';

% Check matric size()

% Compare two spectrogram
sp.height=20; sp.width=40; sp.nbx=1; sp.nby=2;
sp.ledge=3; sp.redge=3.5; sp.tedge=2; sp.bedge=2;
sp.spacex=0.3; sp.spacey=0.7;
%sp.fracy = [0.1 0.3 0.3 0.3];
sp.pos=subplot2(sp);

figure(1)
sizeFont= 24;
set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen

% Creat meash grid
meshNvoie = reshape(1:2,2,1);

for u=1:2
    [sy,sx] = find(meshNvoie==u);
    ax(u)=axes('position',sp.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
    %ax = subplot(2,1,1)
    %subplot(nbSub,nbSub,u)
    %M =  (squeeze(recon(u,:,:)))';
    
    %Plot
    if u==1
    pcolor(t1, f1,Pdb1'); shading flat;
    else
        pcolor(t2, f2,Pdb2'); shading flat;
    end
 
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
    
    title(nameSub{u})
    
    %caxis(spgm.im.clims)
    if u==1
        caxis([30 70])
    else
        caxis([35 74])
    end
    
    
    colormap jet
    %colormap(cmapBW)
    grid on
    
    % Color bar
cbPos = [sp.pos{u}(1)+sp.pos{u}(3)+0.01 sp.pos{u}(2) 0.01 (sp.height/2 - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
ax(u).Position = sp.pos{u};
ylabel(cb,'Power (dB)')
    
end  

%sgtitle([file{i} '  |  ' num2str(floor(timeWin(1))) ' s to ' num2str(ceil(timeWin(end))) ' s  |  File ' num2str(j) '/' num2str(spec.nbIm) ] ,'FontSize', sizeFont+4)



if ~exist('outName')
    outName =['compareSpec.png' ];
end
print('-dpng','-r150',outName)