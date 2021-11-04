% This tool help you to navigate throught file to identify them

% structure
clear all
close all


% Path information : folderIn = wav folder / folderOut = figure output folder
% arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
arrID = 'prc';
outName = 'fileNavigation';
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderIn = ['D:\Bring_Dep_1\' arrID '\']; % Local Mac folder
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];
%folderIn = ['/Volumes/KevDD/BRING/wav/' arrID '/']; % Local Mac folder
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/'];  
if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end


% Loading files and time
% Can specify time = [timeMin, timeMax]
%time = [ min(time), max(time)];

% Load list of file
fileList = getWavName([],folderIn);
fileTime = getFileTime(fileList);
wavInfo = audioinfo([folderIn fileList{1}]);
nbF = length(fileList);

%% Treat wav
i=1;
% Show basic spectro
%[Y,Fs, tstart, dura] = readBring([folderIn fileList{i}], fileTime(i),'duration',wavInfo.TotalSamples/wavInfo.SampleRate,'buffer',0);
[Y,Fs] = audioread([folderIn fileList{i}]);

Y1= Y(:,1);

%%

Nx = length(Y1);
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));

t1 = now;
t = spectrogram(Y1,hamming(nsc),nov,nff);
t2=now;

dt = (t2-t1)

spectrogram(Y1,hamming(nsc),nov,nff);
%% Fourier transform
% Spectro windows
Ns = wavInfo.TotalSamples;
LFFT_spectro = 2048;
LFFT_FV = Ns;
REC = 0.9;
w_pond = kaiser(LFFT_spectro ,0.1102*(180-8.7)); w_pond = w_pond*sqrt(LFFT_spectro/sum(w_pond.^2));
fact_zp =4;

SH =-194;
G = 40;
D = 1;

Y2 = 10^(-SH/20)*10^(-G/20)*D* Y1;

t1=now
[temps, freq, Yf_complexe, Yf_dB] = COMP_STFT_snapshot(Y2,0, Fs, LFFT_spectro, REC, w_pond, fact_zp);
t2=now
dt = t2-t1

Yf = f
%%
figure(1)
pcolor(vec_temps, vec_freq, MAT_t_f_STFT_dB'); shading flat; caxis([30 90]); colorbar
ylim([20 300])
xlabel(' t (s)')
ylabel(' f (Hz)')
title(' Channel 1, dB re. 1µPa2/Hz','interpreter','tex')
set(gca,'FontSize',16)
colormap jet