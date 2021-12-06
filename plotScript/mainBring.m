% This script run beamforming for a any configurations. 
% 
% Config script need to specify path, time, arrayID, spectogram parameter.
%
% Time to process might be specifiy by ptime = [datetime(xxx),datetime(xxx)]
% or by a get function [ploc, ptime ]  = getPingLoc('aav');
%
% Figures to plot can be specify by the argument showFig = [1 3 4..]
% Figure code and description are located in showFigBring.m and
% showGlobalFig.m
%
% last update 03/12/2021 by @kevDuquette


% -------------- Fixed variables --------------------

% Loading file information
[fileList, wavID] = getWavName(ptime, folderIn);
wavInfo = audioinfo([folderIn fileList{1}]);
nbF = length(fileList);

% Array information
[arrLoc, arr] = getArrInfo(arrID);

% Spectro parameters
% Modification of variables into a structure spec. If you get error related
% to those name plase do a ctrl+f and modify new name.
%spec.winSz = 2048;  % LFFT_spectro
%spec.rec = 0.9; % REC
%spec.ovlp = 1-spec.rec;
%spec.wpond = kaiser(spec.winSz ,0.1102*(180-8.7)); 
%spec.wpond = spec.wpond*sqrt(spec.winSz/sum(spec.wpond.^2)); %w_pond
%spec.zp =4; % fact_zp
%spec.fmin = fmin_int;
%spec.fmax = fmax_int;
%spgm.im.fmin = fmin_int;
%spgm.im.fmax = fmax_int;

 % Length of wav in second to load
%duraNs = Ns / 10000;     

% Pcolor plot
%spec.Lmin = 30; % Lmin
%spec.Lmax = 70; % Lmax

% Fix parameter
nbAzimut = 360;      % Number of azimut for calcul
azimut360 = linspace(0,360-360/nbAzimut,nbAzimut);
convPW.SH =-194;      % Parameter for power convertion
convPW.G = 40;
convPW.D = 1;
%duree_film = 60;
c = 1475;           % Sound velocities


% Creating a circle vecteur
%xc =  arr.xh;
%yc = arr.yh;

%% File loop

% Pre-Initialisation
matEnergie = nan(nbF,nbAzimut);
angleA = nan(nbF,nbPk);

for ifile=1:length(fileList)
    disp(['Executing beamforming for file ' num2str(ifile) '/' num2str(length(fileList))])
    close all
    
    % Reading wav and extract values
    [wav.s,wav.fs, time_s, audioInfo] = readBring([folderIn fileList{ifile}], ptime(ifile),'duration',spgm.im.dur,'buffer',buffer);
    %file_wav = fileList{ifile};     % file name alone
    spgm.fs = wav.fs; spgm.ns =  audioInfo.Ns; spgm.im.ns= spgm.ns; spgm.im.dur = audioInfo.dura; 
    spgm = getSpgmWin(spgm);        % Get spectograme windows parameter 
   

    % ------------------ BEAMFORMIGN -----------------------
    wav.db =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;
    %aa = size(wav.db);
    %Nsample = aa(1,1);
    %Nc =aa(1,2);
    %Duree = (Nsample-1)*1/fe;
    %time = (0:1:Nsample-1)*1/fe;
    %t=time; % Need to be find and rename
    
    % Get ponderation matrice
    matPondahf = makePond(arrID,azimut360,spgm.im.ns,spgm.fs,'fmin',spgm.im.fmin,'fmax',spgm.im.fmax);
    
    % Make the beamforming for all direction
%   [Pbf, timebf, freqbf, wavBF] = beamForming(arrID, wav.db , azimut360, spgm,'specmethod','spectro');
    

    %  -------------- Goniométrie ---------------------

    % Set frequence value
    freqAll = (0:1:spgm.ns-1)*spgm.fs/spgm.ns;
    indF = find( (freqAll >= spgm.im.fmin) & (freqAll <= spgm.im.fmax));
    freq = freqAll(indF);
    df = freq(2)-freq(1);
    
    % Create matFFT of each chanel with only freq needed
    matFFT2 = nan(size(wav.db));
    for ii = 1:size(wav.db,2)
        matFFT2(:,ii) = fft(wav.db(:,ii));
    end
    matFFT = matFFT2(indF,:); clear matFFT2;
       
    % Initiate Energie
    Energie = nan(1,length(azimut360));
    
    % Loop for energy over direction
    for u = 1 : length(azimut360)
        matpondhf =squeeze(matPondahf(u,:,:));
        Energie(u) = 0;
        for v = 1 : length(freq)
            vecPond = (matpondhf(:,v)');
            vecFFT  = matFFT(v,:);
            %Energie(u) = Energie(u)+ (abs(sum(vec_pond.*vec_sfft)))^2;
            Energie(u) = Energie(u)+ 1/(spgm.ns*spgm.fs)*(abs(sum(vecPond.*vecFFT)))^2*df;
        end
    end
    
    % Other form of Energie unity
    Energie_NORM = Energie/max(Energie);
    Energie_dB = 10*log10(Energie);
    
    % Find peak of energy
    [pk pkloc] = findpeaks(Energie_NORM,'SortStr','descend','NPeaks', 4 );%'MinPeakHeight',0.3);
    
    % Find the direction of source
    if length(pkloc) < nbPk; pkloc= [pkloc nan(1,nbPk - numel(pkloc)  )] ;end
        
    angleA(ifile, : )  = pkloc;  % Angle of arrival with side lobe
    angleM(1,ifile)  = pkloc(1);  % Max angle of arrival

   
    % Keep the energie value in a matrix
    matEnergie(ifile, :) = Energie;
    
    % Printing figure
    showBring;
end % end loop on file



% Show global figure
%showGlobalFig;


% Saving data
if saveData == true
    if exist('angleR')
        save([folderOut 'dataAngleOfArrival_' outName '.mat'],'angleA','angleM','angleR','ptime','matEnergie')
    else
        save([folderOut 'dataAngleOfArrival_' outName '.mat'],'angleA','angleM','ptime','matEnergie')
    end
end
