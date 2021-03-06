% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette
clear all
close all
%addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
arrID = 'MLB';
% Selected time
ptime= datetime(2021,08,04, 02, 58, 11, 30);
%[ptime, ploc]  = getPingInfo(arrID); % Load ping information
%ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
%ptime = ptime(1);

folderIn = ['/Volumes/BringDD3/Bring_Dep_2/' arrID '/']; % Local Mac folder
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS' ) '_rapTech'];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];
%folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\' outName '\'];
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/prcCircle_Ns14_f150-200hz/'];
%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];

% Figure parameters
showFig = [0];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 4;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = 1.5;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% Creating the output folder
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [50 300];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 80];           % [dB] C limite pcolor
spgm.im.dur = imDur;%'all';         % [s or 'all'] figure duration
%spgm.im.ovlp = 50;                 % [%] image window overlap
spgm.im.figvision = true ;          % [true false] visiblity of figure before saveas jpf file
% spectrogram window parameters
spgm.win.dur = 2048/10000;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 2;                  % [s] spectrogram zero padding
spgm.im.fmin = spgm.im.freqlims(1);spgm.im.fmax = spgm.im.freqlims(2);


%%

% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    %unpackStruct(p);  % Unpack structure to get same workspace as locateBring
else
    % For figure related to locateBring see plotScript/showFigBring.m
    disp(['Executing beamforming for one time at ' datestr(ptime)])
    mainBring;
    %locateBRing;
    
    % Add some figure script located in ../plotScript
    % Show global figure
    %showGlobalFig;
    %showAOACircle; 
    
end


% Saving some pameters

% Loading files and time
% Time must be in datetime format. The file to load will be automatically find
%ptime = datetime(2021,07,15,00,00,00);
%ptime = datetime(2021,07,15,00,00,00):seconds(5):datetime(2021,07,15,01,00,00);
%[ptime, ploc]  = getPingInfo(arrID); % Load ping information
%ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
%ptime = ptime(1);



%% Add some more specified line or figures related to you run
% Enjoy!

% compSpectro(Pdb1,t1,f1,Pdb2,t2,f2,spgm, varargin)
[~,f1,t1,tmp] = spectrogram(wav.pa(:,7),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
Pdb1(:,:) = 10*log10(tmp);

[Pdb22, t2, f2,reconC] = beamForming(arrID, wav.pa , aziCible , spgm,'specmethod','spectro');
Pdb2 = squeeze(Pdb22)';

% Compare two spectrogram
sp.height=20; sp.width=30; sp.nbx=1; sp.nby=2;
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
        pcolor(t1, f1,Pdb1); shading flat;
    else
        pcolor(t2, f2,Pdb2); shading flat;
    end
 
    ylim([spgm.im.fmin spgm.im.fmax])
    set(ax(u),'FontSize',18)
    if sx==1
        %ylabel(' $\mathcal{f}$ (Hz)','interpreter','latex','FontSize',sizeFont)
        %ylabel('$f$ (Hz)','interpreter','latex','FontSize',sizeFont)
        ylabel('f (Hz)')
        %set(gca,'Xtick',[150 160 170 180 190],'xtickLabel',num2cell([150 160 170 180 190]))
    else
        ax(u).YTickLabel = [];
    end
    if sy==2
        %xlabel(' $t$ (s)','interpreter','latex','FontSize',sizeFont)
        xlabel('t (s)')
    else
        ax(u).XTickLabel = [];
    end
    
    if exist('nameSub')
    title(nameSub{u})
    end
    
    %caxis(spgm.im.clims)
    caxis(spgm.im.clims)
    
    colormap jet
    %colormap(cmapBW)
    grid on
    
end  

%sgtitle([file{i} '  |  ' num2str(floor(timeWin(1))) ' s to ' num2str(ceil(timeWin(end))) ' s  |  File ' num2str(j) '/' num2str(spec.nbIm) ] ,'FontSize', sizeFont+4)

% Color bar
cbPos = [sp.pos{end}(1)+sp.pos{end}(3)+0.01 sp.pos{end}(2) 0.01 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
ax(end).Position = sp.pos{end};
ylabel(cb,'PSD (dB re. 1?Pa^2/Hz)','interpreter','tex')
%ylabel(cb,'PSD (dB re. 1$\mu$ Pa${}^2$/Hz)','interpreter','latex','FontSize',sizeFont)

if ~exist('outName')
    outName =['compareSpec.png' ];
end
mkdir2([folderOut])
print('-dpng','-r150',[folderOut '/' outName '_compSpect.png'])




disp('Done !')
