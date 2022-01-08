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

% Fix parameter
nbAzimut = 360;      % Number of azimut for calcul
azimut360 = linspace(0,360-360/nbAzimut,nbAzimut);
convPW.SH =-194;      % Parameter for power convertion
convPW.G = 40;
convPW.D = 1;
%duree_film = 60;
c = 1475;           % Sound velocities

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
   
    wav.pa =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;
    
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
    % OPTIMIZE HERE TO DO NOT CALCULATED ALL FFT
    matFFT2 = nan(size(wav.pa));
    for ii = 1:size(wav.pa,2)
        matFFT2(:,ii) = fft(wav.pa(:,ii));
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
