clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 256;
nCoils = 4;

%% Sensitivity Maps
Mcal = load('data/calibration_4_coils.mat').M;
mask = load('data/mask.mat').mask;
maskCalibration = load('data/maskCalibration.mat').maskCalibration;

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

N = nXPi * nYPi;

noise = sqrt(5e-9)*randn(N,nCoils) + 1i*sqrt(5e-9)*randn(N,nCoils);

MxyPiNoise = MxyPi + noise;

stackKspaceUndersampled = getStackedKspace(MxyPiNoise,nXPi,nYPi,nCoils);

stackKspacePI = zeroFillKspace(stackKspaceUndersampled,R);
stackImagePI = getStackedImage(stackKspacePI);

% visualization of undersampled image
stackImageUndersampled = getStackedImage(stackKspaceUndersampled);
showGrid(abs(stackImageUndersampled),'gray');

% Sense
Gamma = corrcoef(noise);
[imageSense,gmap] = reconstructSense(sensitivityMaps,stackImagePI,Gamma);

% visualize result
plotImage(rescale(abs(imageSense.*maskCalibration)),'gray');
plotImage(rescale(abs(gmap)),'gray');

