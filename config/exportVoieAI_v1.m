%% This script load à wav of 20 hyd and reconstruct 36 directions wav
%
% UPDATES:
% 2021-01-13        Kevin.duquette@dfo-mpo.gc.ca (KD)

clear all
close all
%addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
arrID = 'MLB';

% Time specification
time2Load= datetime(2021,08,04, 00, 53, 00);
%time2Load = datetime(2021,08,04,00,00,00):minutes(5):datetime(2021,08,05,00,00,00);

% Path information | You can also add and use fout = getDirectory('fout')
folderIn = ['~/Documents/MPO/BRing/Data/wav/' arrID '/']; % Local Mac folder
[~, wavi] = getWavName(time2Load, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(time2Load,'yyyymmddTHHMMSS') ];%'MLB_1493_20210804T005254';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/' outName '/' ];

% Figure parameters
showFig = [1 2 3 4 5 6];       % Figure number to print

% Number of direction
nbV =10;

% Load parameter
imDur = 15;
buffer = 0;

% spectrogram image parameters
spgm.im.freqlims = [100 200];       % [Hz] frequency scale boundary limits
spgm.im.clims = [30 70];           % [dB] C limite pcolor
spgm.im.dur = imDur;%'all';                % [s or 'all'] figure duration
spgm.im.ovlp = 50;                   % [%] image window overlap
spgm.im.figvision = true  ;           % [true false] visiblity of figure before saveas jpf file
spgm.im.movm = 100;                  % Movmean parameter for the first panel
spgm.im.FontS = 14;                  % Image font Size
% spectrogram window parameters
spgm.win.dur = 0.18;          % [s] spectrogram window length duration
spgm.win.ovlp = 90;                 % [%] spectrogram window overlap
spgm.win.type = 'kaiser';
spgm.win.opad = 4;                  % [s] spectrogram zero padding

% Other parameter
convPW.SH =-194;
convPW.G = 40;
convPW.D = 1;

% Figure parameter
sp.height=20; sp.width=40; sp.nbx=2; sp.nby=nbV/sp.nbx;
sp.ledge=3; sp.redge=3.5; sp.tedge=2; sp.bedge=2;
sp.spacex=0.3; sp.spacey=0.3;
%sp.fracy = [0.1 0.3 0.3 0.3];
sp.pos=subplot2(sp);

% Automated parameter and variables
if any(showFig); mkdir2(folderOut); end

% List file in the folderIn
fileList = getWavName(time2Load,folderIn);

% Get array information
[aPos arr]  = getArrInfo(arrID);

% Create azimut vector
%azimutV = arr.azimutMax(1):(arr.azimutMax(2)-arr.azimutMax(1))/nbV:arr.azimutMax(2);%(0:nbV-1) * (360 / nbV);
azimutV = arr.azimutMax(1):180/(nbV-1):arr.azimutMax(2);


%% Main 
% Loop on file
for ifile = 1:length(fileList)
    fprintf([datestr(datetime('now')) ' | Load file ' num2str(ifile) '/' num2str(length(fileList)) ' -> ' fileList{ifile} '...']);
    
    wavInfo = audioinfo([folderIn fileList{ifile}]);
    spgm.fs = wavInfo.SampleRate;
    
    % Loop on section of 15 sec on a 5 minute record
    for idir=1%:spgm.im.n
        fprintf('%d.',idir)
    
        % Indice and time related
        ind2Read  = (idir-1)*spgm.im.dur * spgm.fs+1 : idir*spgm.im.dur * spgm.fs;
        timeWin = ind2Read / spgm.fs;
        
        % Reading wav and extract values
        %[wav.s,wav.fs] = audioread([folderIn file{ifile}],[ind2Read(1) ind2Read(end)]);
        [wav.s,wav.fs, wav.t, audioInfo] = readBring([folderIn fileList{ifile}], time2Load,'duration',spgm.im.dur,'buffer',buffer);
        
        % Get spectograme windows parameter
        spgm = getSpgmWin(spgm, audioInfo);        
        
        % Convert wav to pascal
        wav.pa =10^(-convPW.SH/20)*10^(-convPW.G/20)*convPW.D*wav.s;
                
        % Beamforming        
        [Pdb, timeV, freqV, reconWav] = beamForming(arrID, wav.pa , azimutV, spgm);
       
        % Show Figure
        showBeamFormingVisu;

        % Ok Maintenant Xavier tu dois réutiliser la fonction spectrogram
        % avec les parametres spgm pour avec les voie contenue dans
        % reconWav afin de s'assurer que c'est belle et bien un bon siganl
        % wav. Lors que ce sera fait. On pourrait exporter c'est fichier
        % wav dans différent fichier ayant le même non que le fichier
        % inpput mais avec _vAngle(1) -> _vAngle(End) dans el nom du
        % fichier. On va transmettre ceci à Phil par la suite. 
        % Tu pourrais comparer le PSD obtenue avec le PSD de Pdb. Il exite
        % une fonction compSpectro dans plotScript.
        % Amuse toi!
        
    end
    
fprintf('\n')

end


