% This script load the full boat track of Rossboroww aroudn the AAV
% structure
clear all
close all

% Config parameter for locateBring.m
addpath(genpath('../'));

% Need to run data or just open it?
openData = true;

% Path information : folderIn = wav folder / folderOut = figure output folder
% arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
arrID = 'AAV';
outName = 'aavFullLiveBoat_10sec_Ns14_f150-200hz';
folderIn = ['D:\Bring_Dep_1\' arrID '\']; % Local Mac folder
folderOut = ['C:\Users\duquettek\Documents\BRing\results\' arrID '\' outName '\'];
%pingFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircle_7sec_Ns16_f150-200hz/'];
%boatFolder = ['/Users/Administrator/Documents/MPO/BRing/Data/results/AAV/aavCircleBoat_Ns14_f150-200hz/'];


%folderIn = ['Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\' arrID '\'];

% Loading files and time
% Time must be in datetime format. The file to load will be automatically find
g = loadGPXTrack([],arrID,'istart',1,'iend',0,'showFig',2,'printFig',true,'folderout',folderOut);
%b = loadGPXTrack([],arrID);
%ptimeP = getPingInfo(arrID);
%[indexT val ] = getInbetween(ptimeT, [min(ptimeP) max(ptimeP)]);
% Now opening file with a dt specifiy
dt = seconds(10);
ptime = min(g.time):dt:max(g.time);


% Figure parameters
showFig = []       % Figure number to print
saveData = true;   % Save result to .mat
printFig = false;    % Saving figure to a folder
nbPk = 4 ;          % Nomber of side lobe to keep

% Reading parameters
Ns = 2^14;              % Total number of sample
buffer = 0.4;             % Time in second


% Spectro parameter
% Frequence min and max
fmin_int = 150;
fmax_int = 200;

% Get the real angle
%arrLoc = getArrLoc(arrID);
angleR = getRealAngle(arrID, g.lat,g.lon);
g.angleR = angleR;


% If wanted to open alredy run script
if openData == true
    disp('Opening already run data')
    %p = getRunData(pingFolder);
    %b = getRunData(boatFolder);
    bName = dir([folderOut 'data*']);
    load ([folderOut bName.name])
   
else
    locateBRing; % The main loop calculation are locate in this script
    
end


%% Figures

% Show the azimute as fonction of time
showAzigrame;

% Show the matrice of energie with ping location
showMatEngergie;
