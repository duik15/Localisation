% Quick config for only one time analize 
%
% Last update 14/10/21 by @kevDuquette

clear all
close all

% Add all function located in the three to path
addpath(genpath('../'));

% Need to run data or just open it a .mat
openData = false;

% Slected time
ptime= datetime(2021,08,04,00,52,50);

% Path information : folderIn = wav folder / folderOut = figure output folder
arrID = 'MLB';
folderIn = ['\\169.254.116.24\usbshare2-2\Bring_Dep_2\' arrID '\']; % Local Mac folder
[~, wavi] = getWavName(ptime, folderIn);
outName = [arrID '_' wavi.wavID '_' datestr(ptime,'yyyymmddTHHMMSS')];%'MLB_1493_20210804T005254';
folderOut = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\results\' arrID '\' outName '\'];
%folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/prcCircle_Ns14_f150-200hz/'];
%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];

disp(['Executing beamforming for one time at ' datestr(ptime)])

% Loading files and time
% Time must be in datetime format. The file to load will be automatically find
%ptime = datetime(2021,07,15,00,00,00);
%ptime = datetime(2021,07,15,00,00,00):seconds(5):datetime(2021,07,15,01,00,00);
%[ptime, ploc]  = getPingInfo(arrID); % Load ping information
%ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
%ptime = ptime(1);
%ploc = ploc(1,:);


% Figure parameters
showFig = [1 2 3 4 5 6];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep
aziCible =  'max';        % 'max' of value from 0 - 360. The azimut wanted for reconsruction 

% Reading parameters
Ns = 3 * 10000;%2^15 %+ (0.4 * 10000);              % Total number of sample
buffer = 0.5;             % Time in second to add before the ptime

% Spectro parameter
% Frequence min and max
fmin_int = 40;
fmax_int = 100;

% Get the real angle and distance from center
arrLoc = getArrInfo(arrID);
%angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    %unpackStruct(p);  % Unpack structure to get same workspace as locateBring
else
    % For figure related to lcoateBruing see showBring.m
    locateBRing; % The main loop calculation are locate in this script
    
    % Add some figure script located in ../plotScript
    % Show global figure
    %showGlobalFig;
    %showAOACircle; 
    
end


%% Add some more specified line or figures related to you run
% Enjoy!

disp('Done !')
