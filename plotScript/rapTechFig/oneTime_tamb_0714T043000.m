% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette
clear all
close all
%addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
%arrID = 'CLD';
arrID = 'AAV';

% Selected time
%ptime= datetime(2021,07,16, 06, 19 , 03 );
ptime= datetime(2021,07,14, 04, 33, 08);

% Path information | You can also add and use fout = getDirectory('fout')
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
folderIn = ['/Volumes/BringDD3/Bring_Dep_1/' arrID '/']; % Local Mac folder
%folderIn = ['/Volumes/BringDD3/Bring_Dep_1/' arrID '/']; % Local Mac folder
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
%outName = [arrID '_' wavi.wavID '_mom' ];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];

% Figure parameters
showFig = [1]%2 3 4 5 ];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
openData = false; % Need to run data or just open it a .mat
imDur = 15;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = imDur/2;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% Creating the output folder
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [0 500];       % [Hz] frequency scale boundary limits
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
%% 

% Move folder name to add azimut
%if imDur <= 4
movefile(folderOut,[folderOut(1:end-1) '_' num2str(angleM) 'd'])
%else
%    movefile(folderOut,[folderOut(1:end-1) '_' num2str(imDur) 's'])
%end

%% Add some more specified line or figures related to you run
% Enjoy!

i_ch = 1:2:2*nbSP;

for ii=1:nbSP
    [~,freq,time,tmp] = spectrogram(wav.pa(:,i_ch(ii)),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    if ii==1; PSD2 = nan([nbSP size(tmp)]) ; end
    PSD2(ii,:,:) = 10*log10(tmp);
end

%%
% figure 1 : Show 10 chanel
sp.height=30; sp.width=30; sp.nbx=1; sp.nby=5;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.2; sp.spacey=0.2;
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');
%mesh2= mesh2';
nbSP = sp.nbx * sp.nby;
fontS = 12;


fig = figure(1);
movegui(fig,'northwest');
% --------------- Figure --------------------------
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperSize', [sp.width sp.height]);
set(fig, 'PaperPositionMode', 'manual');

for ii=1:nbSP
    [sx,sy] = find(mesh2'==ii);
    % Get the spectogram of first chanel
    
    PSD = squeeze(PSD2(ii,:,:));
    ax(ii)=axes('position',sp.pos{sx,sy},'FontSize',fontS,'Box','on','Layer','top');
    pcolor(time,freq,PSD);
    shading flat;
    colormap(jet)
    %axis tight;
    ylim(spgm.im.freqlims);
    %set(gca,'YScale',spgm.im.freqscale);
    caxis([spgm.im.clims])

    
    % Text and label
    if ii<=5
        text(14.5, spgm.im.freqlims(2) * 0.93,['ch.' num2str(i_ch(ii))], 'backgroundcolor', [1 1 1],'margin',2)
    else
        text(14.5, spgm.im.freqlims(2) * 0.93,['ch.' num2str(i_ch(ii))], 'backgroundcolor', [1 1 1],'margin',2)
    end
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2 + 5,['ch.'])
    %text(timeP(end) *1.02 , spgm.im.freqlims(2)/2-5,[num2str(mainCh(ii))])
    if sx == 1
        yl = ylabel('f (Hz)','FontSize',fontS);
    else
        ax(ii).YTickLabel = [];
    end
    if sy == sp.nby
        xl = xlabel('t (s)','FontSize',fontS) ;
    else
        ax(ii).XTickLabel = [];
    end
    
    set(gca,'layer','top','XMinorTick','on','YMinorTick','on');
    grid on
    box on
    
end



% Color bar
cbPos = [sp.pos{sx,sy}(1)+sp.pos{sx,sy}(3)+0.01 sp.pos{sx,sy}(2) 0.02 (sp.height - sp.bedge - sp.tedge)/sp.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
cb.Position(3) = 0.02;
%ax(end).Position = sp.pos{end};
ylabel(cb,'PSD (dB re. 1µPa^2/Hz)')

% General text and label
%sgtitle([arrID ' | ' datestr(ptime)],'Interpreter','none');

print('-dpng','-r150',[folderOut 'fig1_rapTech_showChanels.png'])



disp('Done !')
