clear all
close all
arrID = 'AAV';
%ptime= datetime(2021,07,15,18,10,58);

[ptime, ploc]  = getPingInfo(arrID); % Load ping information
ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
ptime = ptime(1);

folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
[file, wavi] = getWavName(ptime, folderIn);
outName = ['testDB_' arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS')];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];
imDur = 2;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = 0.5;             % Time in second to add before the ptime
arrLoc = getArrInfo(arrID);
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end
%% Spectrogram parameter
% spectrogram image parameters
spgm.im.freqlims = [20 300];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 100];           % [dB] C limite pcolor
spgm.im.dur = imDur;%'all';         % [s or 'all'] figure duration
spgm.im.figvision = true ;          % [true false] visiblity of figure before saveas jpf file
% spectrogram window parameters
spgm.win.dur = 2048/10000;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 2;                  % [s] spectrogram zero padding
spgm.im.fmin = spgm.im.freqlims(1);spgm.im.fmax = spgm.im.freqlims(2);

% Parameter for power convertion 
convPW.SH =-194;  convPW.G = 40; convPW.D = 1;

%% Opening thwe file
% Reading wav and extract values
[wav.s,wav.fs, time_s, audioInfo] = readBring([folderIn file{1}], ptime(1),'duration',spgm.im.dur,'buffer',buffer);
wav.db =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;

% Getting some information
spgm.fs = wav.fs; spgm.ns =  audioInfo.Ns; spgm.im.ns= spgm.ns; spgm.im.dur = audioInfo.dura;
spgm = getSpgmWin(spgm);        % Get spectograme windows parameter


mat2send = wav.db;

%% Version Spectogram
% Get the spectogram of first chanel
    [~,freqS,timeS,tmp] = spectrogram(mat2send(:,1),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    PdbS = 10*log10(tmp);

%% Version Cedric
%function [vec_temps, vec_freq, MAT_t_f_STFT_complexe, MAT_t_f_STFT_dB] = COMP_STFT_snapshot(vec_s, t0, fe, Lfft, REC, w_pond, fact_zp)
%[timeV, freqFV, ~ , matFFTV_dB] = COMP_STFT_snapshot(recon(:,u),0, spgm.fs, spgm.win.ns, spgm.win.ovlp/100, spgm.win.val, spgm.win.opad);
% Debug fonction argument
vec_s = mat2send(:,1);
t0 =1; fe = wav.fs; Lfft = spgm.win.ns; REC = spgm.win.ovlp/100; w_pond = spgm.win.val; fact_zp = spgm.win.opad;

%
Nd1_rel=Lfft/2+1; % indice de début
Nf1_rel=length(vec_s)-Lfft/2;% indice de fin
Pas_c=floor(Lfft*(1-REC)); % pas entre deux fenetres FFT pour assurer un taux de recouvrement == REC
vecteur_indice_centre=Nd1_rel:Pas_c:Nf1_rel;
freq_dsp=(0:fact_zp*Lfft-1)*fe/(fact_zp*Lfft);
indi_int = find( (freq_dsp>=0) & (freq_dsp<fe/2));
vec_freq =freq_dsp(indi_int);
vec_temps = vecteur_indice_centre*1/fe + t0;
%h = waitbar(0,'CALCUL DU SPECTROGRAMME IN PROGRESS');
for u = 1 : length(vecteur_indice_centre) % on parcourt les snapshots
    %waitbar(u/length(vecteur_indice_centre))
    ss_loc  = vec_s(vecteur_indice_centre(u)-Lfft/2+1:vecteur_indice_centre(u)+Lfft/2); % on vient sélectionner le snapshot interessant
    a = size(ss_loc);
    aa = size(w_pond);
    if a(1,1)~=aa(1,1)
        w_pond=w_pond';
    end
    fft_loc = fft(ss_loc.*w_pond,fact_zp*Lfft);
    fft_loc_int = fft_loc(indi_int);
    
    if iscolumn(fft_loc_int)
        MAT_t_f_STFT_complexe(u,:)=fft_loc_int';
    else
        MAT_t_f_STFT_complexe(u,:)=fft_loc_int;
    end
end
%close(h)
timeC = vec_temps;
freqC = vec_freq;
PdbC = 10*log10(1/(Lfft*fe)*abs(MAT_t_f_STFT_complexe).^2);

%% 
    
    % Get the spectogram of first chanel
    [~,freq1,time1,tmp] = spectrogram(wav.db(:,1),spgm.win.val,spgm.win.novlp,spgm.win.nfft,spgm.fs);
    Pdb1 = 10*log10(tmp);
    
    figure(2)
    pcolor(timeS, freqS, PdbS); shading flat;
    ylim([20 300])
    xlabel(' t (s)')
    ylabel(' f (Hz)')
    title(' Channel 1, dB re. 1µPa2/Hz','interpreter','tex')
    set(gca,'FontSize',16)
    colormap jet
    caxis(spgm.im.clims)
    cb = colorbar;
    ylabel(cb,'Power (dB)')


%% Compare spectro
compSpectro(PdbS,timeS,freqS,PdbC',timeC,freqC,spgm)