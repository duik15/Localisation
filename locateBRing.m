% This script aim to locate the structure from our upcall emission
clear all

% Import WAV
folderIn = '../wavBRING/';%'/Volumes/MartinUSB/'; 'Z:\DATA\missions\2021-07-27_IML_2021-016_BRings\wav\';   
fileList = {'2591LF_20210714T222000_0557_0.wav'};
wavInfo = audioinfo([folderIn fileList{1}]);
[Y,Fs] = audioread([folderIn fileList{1}]);



%plot(Y)
