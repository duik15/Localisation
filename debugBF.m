clc
clear all
close all
tic
%% Variables
% enviroment variables
arrID = 'MLB';
%outName = 'aavLiveBoat_onlyTrackTime_Ns14_f150-200hz';
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
folderOut = [folderIn '/' 'beamFormingAll/'];
outName = '';
%folderIn = ['E:\Bring_Dep_1\' arrID '\'];
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\spectro\']; %spectrogram'];

saveflag = false;                   % [true or false]

% Selecting file wanted
time2Load = datetime(2021,08,04,00,55,00);
%time2Load = getPingInfo('AAv'); time2Load = time2Load(1);
%time2Load = datetime(2021,07,14,00,00,00):minutes(5):datetime(2021,07,14,02,00,00);

% Number of direction
nbV =10;


% Frequence min and max
fmin = 100; fmax =200;

% List file in the folderIn
file = getWavName(time2Load,folderIn);

% Get array information
[aPos arr]  = getArrInfo(arrID);

% Create azimut vector
azimutV = arr.azimutMax(1):(arr.azimutMax(2)-arr.azimutMax(1))/nbV:arr.azimutMax(2);%(0:nbV-1) * (360 / nbV);

if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%%

        
        
%%
% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.freqscale = 'linear';       % ['linear' or 'log'] frequency scale type
spgm.im.clims = [-120 0];           % [dB] C limite pcolor
spgm.im.lbwh = 'auto';              % [left, bottom, width, height or 'auto'] image size and location [800 500 188 179];
spgm.im.format = 'png';             % image output format 'png', 'jpg'
spgm.im.dur = 15;%'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 0;                   % [%] image window overlap
spgm.im.link_axes = false;           % [true false] use fonction linkaxes
spgm.im.figvision = true  ;           % [true false] visiblity of figure before saveas jpf file
spgm.im.constrast = [1 1.2];        % Constrat the clims by [x x] between 0 and 1 caxis = caxis * [cont1 cont2]
spgm.im.movm = 100;                  % Movmean parameter for the first panel
spgm.im.FontS = 14;                  % Image font Size
% spectrogram window parameters
spgm.win.dur = 0.18;                  % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 4;                  % [s] spectrogram zero padding

spec.winSz = 0.18 * 10000;  % LFFT_spectro
spec.rec = 0.9; % REC
sepc.ovlp = 1-spec.rec;
spec.wpond = kaiser(spec.winSz ,0.1102*(180-8.7));
spec.wpond = spec.wpond*sqrt(spec.winSz/sum(spec.wpond.^2)); %w_pond
spec.zp =4; % fact_zp
spec.fmin = fmin;
spec.fmax = fmax;
spec.Lmin = 30; % Lmin
spec.Lmax = 70; % Lmax
spec.SH =-194;
spec.G = 40;
spec.D = 1;
%%
i=1
j=1
acinfo = audioinfo([folderIn file{1}]);

       ind2Read  = (j-1)*spgm.im.dur * acinfo.SampleRate+1 : j*spgm.im.dur * acinfo.SampleRate;
        timeWin = ind2Read / acinfo.SampleRate;
        
acinfo = audioinfo([folderIn file{i}]);
        [wav.s,wav.fs] = audioread([folderIn file{i}],[ind2Read(1) ind2Read(end)]);
        wav.ns = size(wav.s,1);
        wav.ch = size(wav.s,2);
        % Power conversion
        wav.sdb = 10^(-spec.SH/20)*10^(-spec.G/20)*spec.D* wav.s;
        
 %%
        
acf = wav;

acinfo.SampleRate = acf.fs
spgm.im.ns = fix(spgm.im.dur*acinfo.SampleRate);
spgm.im.n = round(acinfo.TotalSamples/spgm.im.ns);

%Nx = length(x);NS
%nsc = floor(Nx/4.5);
%nov = floor(nsc/2);
%nff = max(256,2^nextpow2(nsc));
%t = spectrogram(x,hamming(nsc),nov,nff);

% define spectrogram window
win.ns = fix(spgm.win.dur*acf.fs);
win.wpond = kaiser(win.ns ,0.1102*(180-8.7));
win.val = win.wpond*sqrt(win.ns/sum(win.wpond.^2)); %w_pond
win.novlp = fix(spgm.win.ovlp/100*win.ns);
win.nfft = win.ns*spgm.win.opad;
%win.val = %hanning(win.ns);
%win.nfft = fix((spgm.win.dur+ spgm.win.opad)*acf.fs);  %max(256,2^nextpow2(win.ns*10))


% need to remove
spec.Ns = spgm.im.ns
Ns = spec.Ns
spec.fmin = spgm.im.freqlims(1);
spec.fmax = spgm.im.freqlims(2);
spec.Fs = wav.fs;
spec.nbIm = 3333;
%
% Get the ponderation matrice
matPondahf = makePond(arrID,azimutV,spec.Ns,spec.Fs,'fmin',spec.fmin,'fmax',spec.fmax);


%%
matWav = wav.sdb(:,10);
% Create the frequence vecteur relate to FFT
freq= (0:1:Ns-1)*spec.Fs/Ns;
indF = find( (freq >= spec.fmin)&(freq <= spec.fmax));
freqF = freq(indF);

% Create matFFT with only freq needed
matFFT2 = nan(size(matWav));
for ii = 1:size(matWav,2)
    matFFT2(:,ii) = fft(matWav(:,ii));
end
matFFT = matFFT2(indF,:); clear matFFT2;

%% Initiate
clear reconFFT reconFFTC
recon=zeros(spec.Ns,nbV);

u=1
    matFFTPond = zeros(1,floor(Ns/2)+1);
    % Calcul the fft for the direction and add those freq values to the null matrix
    matFFTPondSum = ( sum(matFFT .* squeeze(matPondahf(u,:,:))',2 ))';
    
    % Add the conj
    matFFTPond(indF) = matFFTPondSum;
    matFFTPond2 = [ matFFTPond conj(fliplr(matFFTPond(2:length(matFFTPond)-1))) ];
    %matFFTPondV(:,u) =  matFFTPond2';
    % Reconstruct WAV % This matrice could be remove to save memory and
    % time
    recon(:,u) = ifft(matFFTPond2);
    
    % Calcul of spectrogrammes associated with BF
    [~,freqV,timeV,mat1] = spectrogram(recon(:,1),win.val,win.novlp,win.nfft,spec.Fs); ok1=1
    matFFTV_dB = 10*log10(mat1');
    %[~,freqV,timeV,mat1] = spectrogram(matWav(:,1),win.val,win.novlp,win.nfft,spec.Fs); ok1=1
    %[timeV2, freqV2, ~ , mat2] = COMP_STFT_snapshot(matWav(:,1),0, spec.Fs, spec.winSz, spec.rec, spec.wpond, spec.zp);ok2=1
    %freqV=freqV2; timeV=timeV2
    %matFFTV_dB = mat2';
    [timeVC, freqVC, ~ , matFFTV_dBC] = COMP_STFT_snapshot(recon(:,u),0, spec.Fs, spec.winSz, spec.rec, spec.wpond, spec.zp);
    
    %[~,freqV,timeV,matFFTV_dB] = spectrogram(recon(:,u),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spec.Fs); 
    %[timeV, freqV, ~ , matFFTV_dB] = COMP_STFT_snapshot(recon(:,u),0, spec.Fs, spec.winSz, spec.rec, spec.wpond, spec.zp);
    indFr = find( ((freqV>= spec.fmin)&(freqV <= spec.fmax)));
    indFrC = find( ((freqVC>= spec.fmin)&(freqVC <= spec.fmax)));

    %reconFFT(u,:,:)=matFFTV_dB(indFr,:)';
    reconFFT(u,:,:)=matFFTV_dB(:,indFr);
    reconFFTC(u,:,:)=matFFTV_dBC(:,indFrC);
    
    freqVAll = freqV;
    freqVAllC = freqVC;
    freqV = freqVAll(indFr);
    freqVC = freqVAllC(indFrC);
    disp('Done!')
    %
     %%
     disp('new')
     st = size(timeVC)
     sf = size(freqVC)
     sm = size((squeeze(reconFFTC(u,:,:)))')
     
     %%
     sp2.height=20; sp2.width=40; sp2.nbx=1; sp2.nby=2;
sp2.ledge=3; sp2.redge=3.5; sp2.tedge=2; sp2.bedge=2;
sp2.sp2acex=0.3; sp2.sp2acey=0.3;
%sp2.fracy = [0.1 0.3 0.3 0.3];
sp2.pos=subplot2(sp2);

figure(1)
sizeFont= 24;
set(gcf, 'PaperPosition', [0 0 20 10]);    % can be bigger than screen

% Creat meash grid
meshNvoie = reshape(1:2,2,1);


for u=1:2
    [sy,sx] = find(meshNvoie==u);
    ax(u)=axes('position',sp2.pos{sx,sy},'XGrid','off','XMinorGrid','off','FontSize',14,'Box','on','Layer','top');
    %ax = subplot(2,1,1)
    %subplot(nbSub,nbSub,u)
    %M =  (squeeze(recon(u,:,:)))';
    
    %Plot
    if u==1
    pcolor(timeVC + timeWin(1), freqVC,(squeeze(reconFFTC(1,:,:)))'); shading flat;
    else
        pcolor(timeV + timeWin(1), freqV,(squeeze(reconFFT(1,:,:)))'); shading flat;
       %pcolor(timeV2 + timeWin(1), freqV2,(squeeze(reconFFT2(u,:,:)))'); shading flat; 
    end
    ylim([spec.fmin spec.fmax])
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
    
    if u <= nbV/2
        tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        tt.Position = [ timeWin(1)-0.12*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1) ];
    else
        tt = title([ num2str(azimutV(u)) '$^{\circ}$'],'FontSize',sizeFont,'interpreter','latex');
        tt.Position = [ timeWin(end)+0.15*(timeWin(end)-timeWin(1))  (round(freqV(end)-freqV(1)))/2+freqV(1)-5 ];
    end
    
    caxis([spec.Lmin spec.Lmax])
    colormap jet
    %colormap(cmapBW)
    grid on
    %%
end  

sgtitle([file{i} '  |  ' num2str(floor(timeWin(1))) ' s to ' num2str(ceil(timeWin(end))) ' s  |  File ' num2str(j) '/' num2str(spec.nbIm) ] ,'FontSize', sizeFont+4)

% Color bar
cbPos = [sp2.pos{end}(1)+sp2.pos{end}(3)+0.01 sp2.pos{end}(2) 0.01 (sp2.height - sp2.bedge - sp2.tedge)/sp2.height];%sp.pos{end}(4)*(sp.nby-1)+sp.pos{1}(4)+sp.spacey ];
cb= colorbar;%('Position',cbPos);
cb.Position = cbPos;
ax(end).Position = sp2.pos{end};
ylabel(cb,'Power (dB)')

print('-dpng','-r150',[folderOut 'debug_BF_' file{i}(1:end-4) '_s' num2str(floor(timeWin(1))) '_e' num2str(ceil(timeWin(end))) '.png' ])