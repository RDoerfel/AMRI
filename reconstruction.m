clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 256;
nCoils = 4;

%% Sensitivity Maps
Mcal = load('data/calibration_8_coils.mat').M;
mask = load('data/mask.mat').mask;
maskCalibration = load('data/maskCalibration.mat').maskCalibration;

MxyCal = getMxy(Mcal);

stackKspaceCal = getStackedKspace(MxyCal,nXCal,nYCal,nCoils);

stackImageCal = getStackedImage(stackKspaceCal);
showGrid(abs(stackImageCal).*maskCalibration,'gray');

sensitivityMaps = computeCoilSensitivities(stackImageCal);
sensitivityMapsSmothed = smootheSensitityMaps(sensitivityMaps,1);
showGrid(abs(sensitivityMapsSmothed),'gray');

%% PI
nXPi = 256;
nYPi = 128;
Mpi = load('data/r2_8coils.mat').M;

MxyPi = getMxy(Mpi);

R = nYCal/nYPi;
stackKspaceUndersampled = getStackedKspace(MxyPi,nXPi,nYPi,nCoils);

stackKspacePI = zeroFillKspace(stackKspaceUndersampled,R);
stackImagePI = getStackedImage(stackKspacePI);

% visualization of undersampled image
stackImageUndersampled = getStackedImage(stackKspaceUndersampled);
showGrid(abs(stackImageUndersampled),'gray');

% Sense
Gamma = eye(nCoils,nCoils);
[imageSense,gmap] = reconstructSense(sensitivityMapsSmothed,stackImagePI,Gamma,R);

% visualize result
plotImage(rescale(abs(imageSense.*maskCalibration)),'gray');
plotImage(medfilt2(abs(gmap.*maskCalibration)),'gray');
colorbar;
