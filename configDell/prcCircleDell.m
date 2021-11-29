clear all
close all

% Config parameter for locateBring.m
addpath(genpath('../'));


% Need to run data or just open it?
openData = false;


% Path information : folderIn = wav folder / folderOut = figure output folder
% arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
arrID = 'PRC';
outName = 'charlotte_prcCircle_Ns14_f150-200hz';
folderIn = ['F:\Bring_Dep_2\' arrID '\']; % Local Mac folder
folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];    %folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];
% Loading files and time
% Time must be in datetime format. The file to load will be automatically find
[ptime, ploc]  = getPingInfo(arrID); % Load ping information
ptime = ptime ;%+ seconds(7);          % Add an offset to be center the 5 upcalls

% Figure parameters
showFig = [1 2 3 4 7]       % Figure number to print
saveData = true;   % Save result to .mat
printFig = true;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^14;              % Total number of sample
buffer = 0.4;             % Time in second


% Spectro parameter
% Frequence min and max
fmin_int = 150;
fmax_int = 200;


% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    
else
    locateBRing;
    
    showGlobalFig;;
    
end