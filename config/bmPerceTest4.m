% This script make a triangularization of two beam to find back to position

% structure
clear all
close all

% Config parameter for locateBring.m
addpath(genpath('../'));

% Path information : folderIn = wav folder / folderOut = figure output folder
% arrID = AAV / CLD / MLB / PRC || outName = name of the output folder and figures
arrID = {'mlb', 'prc'};
outName = 'percLiveBoat_beamcrossing_60sec_Ns14_f150-200hz';
folderOut = ['/Users/Administrator/Documents/MPO/BRing/Data/results/perce/' outName '/'];  

inName = {'mlbFullLiveBoat_late_60sec_Ns14_f150-200hz/','prcFullLiveBoat_fitMLB_60sec_Ns14_f150-200hz/'};

% Load information
for i=1:2
folderIn = ['/Users/Administrator/Documents/MPO/BRing/Data/results/' arrID{i} '/' inName{i}]; % Local Mac folder
dataName = dir([folderIn 'data*']);
p{i} = load ([folderIn dataName.name])
arrLoc(i,:) = getArrLoc(arrID{i});
end

% Angle for testing
ia =4;
aa = [p{1}.angleM(ia) p{2}.angleM(ia)];

%%
pos1 = arrLoc(1,:);
pos2 = arrLoc(2,:);
a1 = aa(1);
a2 = aa(2);