%% WAV to fig beam formig spectrogramme
%
% DESCRIPTION:
% Load audio file, make create a beamforming and convert fig in jpg file
% Use IML-PAM-Toolbox\wav_tools\wav_spectrogram_image.m to make this script developed by Florian A and Clement J
%
%
% UPDATES:
% 2021-11-20        Kevin.duquette@dfo-mpo.gc.ca (KD)
%


clc
clear all
close all
tic
%mem0 = memory;
%% Variables
% enviroment variables
arrID = 'PRC';
%outName = 'aavLiveBoat_onlyTrackTime_Ns14_f150-200hz';
%folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
%folderOut = [folderIn '/' 'beamFormingAll/'];
outName = '';
folderIn = ['\\169.254.47.215\usbshare1-2\Bring_Dep_2\' arrID '\']
folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\beamForming\']; %spectrogram'];

saveflag = false;                   % [true or false]

% Selecting file wanted
%time2Load = datetime(2021,07,15,18,05,00);
time2Load = datetime(2021,08,04,00,00,00):minutes(5):datetime(2021,08,05,00,00,00);

% Number of direction
nbV =10;


% Frequence min and max
fmin = 150; fmax =200;

% List file in the folderIn
file = getWavName(time2Load,folderIn);

% Get array information
[aPos arr]  = getArrInfo(arrID);

% Create azimut vector
azimutV = arr.azimutMax(1):(arr.azimutMax(2)-arr.azimutMax(1))/nbV:arr.azimutMax(2);%(0:nbV-1) * (360 / nbV);

if ~isfolder(folderOut); disp(['Creating output folder: ' folderOut]); mkdir(folderOut); end

% Figure parameter
sp.height=20; sp.width=40; sp.nbx=2; sp.nby=nbV/sp.nbx;
sp.ledge=3; sp.redge=3.5; sp.tedge=2; sp.bedge=2;
sp.spacex=0.3; sp.spacey=0.3;
%sp.fracy = [0.1 0.3 0.3 0.3];
sp.pos=subplot2(sp);

% Spectrogramme and beamforming information
% Spectro parameters
% Modification of variables into a structure spec. If you get error related
% to those name plase do a ctrl+f and modify new name.
spec.winSz = 2048;  % LFFT_spectro
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
spec.dur = 15;

spec.nbIm = 300 / spec.dur;
nbIm =  spec.nbIm;
% Other variables needed
%duraNs = Ns / 10000;
c = 1475;    % Sound velocity


%% Open the first

acinfo = audioinfo([folderIn file{1}]);
winL = spec.dur * acinfo.SampleRate;
[wav.s,wav.fs] = audioread([folderIn file{1}],[1 winL]);
spec.Fs = wav.fs;
spec.Ns = size(wav.s,1);

%acf.s = acf.s(:,mainCh);
if strcmp(spec.dur,'all')
    spec.dur = acinfo.Duration;
end


%% process
tic
% Get the ponderation matrice
matPondahf = makePond(arrID,azimutV,spec.Ns,spec.Fs,'fmin',spec.fmin,'fmax',spec.fmax);

% Loop on file
for i = 1:length(file)
    % Read file
    %disp([datestr(datetime('now')) ' | Load file ' num2str(i) '/' num2str(length(file)) ' -> ' file{i}]);
    fprintf([datestr(datetime('now')) ' | Load file ' num2str(i) '/' num2str(length(file)) ' -> ' file{i} '...']);
    % Loop on windows
    for j=1:nbIm
        fprintf('%d.',j)
    
        % Indice and time related
        ind2Read  = (j-1)*spec.dur * acinfo.SampleRate+1 : j*spec.dur * acinfo.SampleRate;
        timeWin = ind2Read / acinfo.SampleRate;
        
        % Load audio
        acinfo = audioinfo([folderIn file{i}]);
        [wav.s,wav.fs] = audioread([folderIn file{i}],[ind2Read(1) ind2Read(end)]);
        wav.ns = size(wav.s,1);
        wav.ch = size(wav.s,2);
        % Power conversion
        wav.sdb = 10^(-spec.SH/20)*10^(-spec.G/20)*spec.D* wav.s;
        
        % Beamforming
        [reconFFT, timeV, freqV] = beamForming(matPondahf, wav.sdb , azimutV, spec);
        
        % Show Figure
        showBeamFormingVisu;
        
        % Save time statistic
        tictocPrint(i,j) = toc;
        close all
        
    end
    
fprintf('\n')
tictocFile(i) = toc;
end
tictocScript = toc;

% Statistic
totalIm = length(file) * nbIm;
diffPrint=diff(tictocPrint);
diffFile=diff(tictocFile);
nowName = datestr(datetime('now'),'yyyymmddTHHMMSS');
% Print stats
save([folderOut 'tictoc_' nowName '.mat'], 'file', 'nbIm','tictocPrint','tictocFile','tictocScript','totalIm','diffPrint','diffFile','nowName')

