% Defautl config for using the beamforming for BRING
%
% Last update 14/10/21 by @kevDuquette

clear all
close all

% Add all function located in the three to path
addpath(genpath('../'));

% Need to run data or just open it a .mat
openData = false;

% Path information : folderIn = wav folder / folderOut = figure output folder
arrID = 'AAV';
outName = 'aav_Ns16_Voc1-5_fm150-200Hz';
folderIn = ['D:\Bring_Dep_1\' arrID '\']; % Local Mac folder
folderOut = ['C:\Users\duquettek\Documents\BRingsOne\results\' arrID '\' outName '\'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID '/prcCircle_Ns14_f150-200hz/'];
%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];

% Loading files and time
% Time must be in datetime format. The file to load will be automatically find
%ptime = datetime(2021,07,15,00,00,00);
%ptime = datetime(2021,07,15,00,00,00):seconds(5):datetime(2021,07,15,01,00,00);
[ptime, ploc]  = getPingInfo(arrID); % Load ping information
ptime = ptime + seconds(7);          % Add an offset to be center the 5 upcalls
%ptime = ptime(1);
%ploc= ploc(1,:);

% Figure parameters
showFig = [5];       % Figure number to print
saveData = true;   % Save result to .mat
printFig = false;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^16;              % Total number of sample
buffer = 0.4;             % Time in second to add before the ptime

% Spectro parameter
% Frequence min and max
fmin_int = 150;
fmax_int = 200;

% Get the real angle and distance from center
arrLoc = getArrLoc(arrID);
angleR = getRealAngle(arrID, ploc(:,1),ploc(:,2));

% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    p = getRunData(pingFolder);
    %unpackStruct(p);  % Unpack structure to get same workspace as locateBring
else
    locateBRing; % The main loop calculation are locate in this script
    
    % Add some figure script located in ../plotScript
    
    %showAOACircle; 
    
end


%% Add some more specified line or figures related to you run
% Enjoy!

