% Defautl config for using the beamforming for BRING
% Mainly you need to specified the:
% 
% arrID = AAV / CLD / MLB / PR 
% outName = name of the output folder and figures
% folderIn = where you wav are
% folderOut = where you want to figure and data to be saved
% ptime = time to locate -> can aslo be used with getPingLoc if ping...
% information are store in the function with and ID related to the data
%
% If openData == False the code will run the locateBring script and saved
% result. If openData = False, the code will load the data matching the
% outName folder in the [folderOut outName] or pingPath that you specified.
% This save calcul time when working many time on the same dataset.
% 
% Other parameter can be set in this script.
%
% Data in you folder should be  
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
outName = 'defaultConfig';
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
ptime = ptime(1);
ploc = ploc(1,:);

% Figure parameters
showFig = [3];       % Figure number to print
saveData = false;   % Save result to .mat
printFig = false;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^14;              % Total number of sample
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

