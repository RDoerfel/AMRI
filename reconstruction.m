clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 128;
nCoils = 4;

%% Sensitivity Maps
Mcal = load('data/calibration_4_coils.mat').M;
load('data/mask.mat');
load('data/maskCalibration.mat');

MxyCal = getMxy(Mcal);

stackKspaceCal = getStackedKspace(MxyCal,nXCal,nYCal,nCoils);
stackImageCal = getStackedImage(stackKspaceCal);

sensitivityMaps = computeCoilSensitivities(stackImageCal);

showGrid(abs(sensitivityMaps).*maskCalibration,'gray');
showGrid(abs(stackImageCal).*maskCalibration,'gray');

%% PI
nXPi = 256;
nYPi = 128;
Mpi = load('data/r2_4_coils.mat').M;

MxyPi = getMxy(Mpi);

R = nYCal/nYPi;

stackKspaceUndersampled = getStackedKspace(MxyPi,nXPi,nYPi,nCoils);
stackKspacePI = zeroFillKspace(stackKspaceUndersampled,R);
stackImagePI = getStackedImage(stackKspacePI);

% visualization of undersampled image
stackImageUndersampled = getStackedImage(stackKspaceUndersampled);
showGrid(abs(stackImageUndersampled),'gray');

%% Sense
[imageSense,gmap] = reconstructSense(sensitivityMaps,stackImagePI);
plotImage(rescale(abs(imageSense.*maskCalibration)),'gray');
plotImage(rescale(abs(gmap)),'gray');

