
clear all
close all


% Path information : folderIn = wav folder / folderOut = figure output folder
arrID = 'MLB';

% Selected time
time= datetime(2021,08,04,00,52,46);

%folderIn = ['F:\Bring_Dep_1\' arrID '\']; % Local Mac folder
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
%folderIn = ['\\169.254.116.24\usbshare2-2\Bring_Dep_2\' arrID '\']
%folderIn = ['C:\Users\duquettek\Documents\BRing\Data\wav\' arrID '\'];
[~, wavi] = getWavName(time, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(time,'yyyymmddTHHMMSS') '_rapTech15s'];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];
%folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\' outName '\caxis3080\'];
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/prcCircle_Ns14_f150-200hz/'];
%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];

% Reading parameters
showFig =2 ;
imDur = 15;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = 0;             % Time in second to add before the ptime

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% Creating the output folder
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%%

% Read
[file, wavID] = getWavName(time, folderIn);
file = file{1};
fileTime = getFileTime(file);

[wav.s,wav.fs, time_s, audioInfo] = readBring([folderIn file], time,'duration',imDur);


convPW.SH =-194;      % Parameter for power convertion
convPW.G = 40;
convPW.D = 1;
wav.pa =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;


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
spgm.fs = wav.fs; spgm.ns =  audioInfo.Ns; spgm.im.ns= spgm.ns; spgm.im.dur = audioInfo.dura; 
spgm = getSpgmWin(spgm);        % Get spectograme windows parameter     


%%
if any(showFig ==1)
sp.height=30; sp.width=30; sp.nbx=1; sp.nby=2;
sp.ledge=4; sp.redge=6; sp.tedge=3; sp.bedge=3;
sp.spacex=0.5; sp.spacey=0.3;
sp.fracy=[0.1 0.9];
[sp.pos, mesh2]=subplot2(sp,'meshorder','x');


%i_ch = 1;
 for i_ch =1 : size(wav.pa,2) 
    close all
    disp(['Printing with amp ' num2str(i_ch) '/' num2str(size(wav.pa,2) ) ])
     
    [~,freq1,time1,tmp] = spectrogram(wav.pa(:,i_ch),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    Pdb1 = 10*log10(tmp);
    
    figure(1)
    ax(1)=axes('position',sp.pos{1,1},'Box','on','Layer','top');
    plot(time_s,wav.s(:,i_ch),'k')
    ylabel('Amp.')
    ax(1).XTickLabel = [];
    %title(['Azimut ' num2str(azimut360(indAziCible)) '°'])
    
    ax(2)=axes('position',sp.pos{1,2},'Box','on','Layer','top');
    
    pcolor(time1, freq1,Pdb1); shading flat;
    ylim([spgm.im.freqlims])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    caxis(spgm.im.clims)
    set(ax(2),'Box','on','Layer','top');
    %caxis([Lmin Lmax])
    colormap jet
    
    % Coordinate axis
    ax(1).XLim = ax(2).XLim;
    
    % Colorbar
    tmppos = ax(2).Position;
    cb = colorbar;
    ax(2).Position = tmppos;
    cb.Position(3) = 0.01;
    ylabel(cb,'PSD (dB re. 1µPa2/Hz)')
    
    print('-dpng','-r150',[folderOut 'showWavChanelsAmp_' outName '_ch' num2str(i_ch) '.png'])
 end    
end
%% No amp
if any(showFig ==2)
clear ax
%i_ch = 1;
 for i_ch = 1 : size(wav.pa,2)   
     close all
     disp(['Printing no amp ' num2str(i_ch) '/' num2str(size(wav.pa,2) ) ])
     
    [~,freq1,time1,tmp] = spectrogram(wav.pa(:,i_ch),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    Pdb1 = 10*log10(tmp);
    
    figure(1)
    pcolor(time1, freq1,Pdb1); shading flat;
    ylim([spgm.im.freqlims])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    caxis(spgm.im.clims)
    %caxis([30 90])
    ax = gca;
    set(ax,'Box','on','Layer','top');
    %caxis([Lmin Lmax])
    colormap jet
    
    set(gca,'XAxisLocation','top')
    
    % Colorbar
    ax(1).Position(1) = 0.1;
    tmppos = ax(1).Position;
    cb = colorbar;
    ax(1).Position = tmppos;
    cb.Position(3) = 0.02;
    ylabel(cb,'PSD (dB re. 1µPa^2/Hz)')
    
    print('-dpng','-r150',[folderOut 'showWavChanels_' outName '_ch' num2str(i_ch) '.png'])
 end    

end
%%
%{
% Get the spectogram of first chanel
    [~,freq1,time1,tmp] = spectrogram(wav.pa(:,1),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    Pdb1 = 10*log10(tmp);
    
    figure(1)
    pcolor(time1, freq1, Pdb1); shading flat;
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1µPa2/Hz','interpreter','tex')
    set(gca,'FontSize',16)
    colormap jet
    caxis(spgm.im.clims)
    cb = colorbar;
    ylabel(cb,'Power (dB)')
%}