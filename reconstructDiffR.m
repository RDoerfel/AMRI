clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 256;
nCoils = 8;

M = load('data/calibration_8_coils.mat').M;
mask = load('data/maskCalibration.mat').maskCalibration;

Mxy = getMxy(M);
stackKspace = getStackedKspace(Mxy,nXCal,nYCal,nCoils);
stackImageCal = getStackedImage(stackKspace);
sensitivityMaps = computeCoilSensitivities(stackImageCal);

%% Different Rs 
MR2 = load('data/r2_8_coils.mat').M;
nXCal = 256;
nYCal = 128;

MxyR2 = getMxy(MR2);
kspaceR2 = getStackedKspace(MxyR2,nXCal,nYCal,nCoils);
imgR2 = getStackedImage(kspaceR2);

R = 2;
Gamma = eye(nCoils,nCoils);
[imgSenseR2,gMapR2] = reconstructSense(sensitivityMaps,imgR2,Gamma,R);
mean(abs(gMapR2(mask)))

% visualize result
plotImage(rescale(abs(imgSenseR2.*mask)),'gray');
plotImage(medfilt2(abs(gMapR2.*mask)),'gray');
colorbar;

% R3
MR3 = load('data/r3_8_coils.mat').M;
nXCal = 256;
nYCal = 128;

MxyR3 = getMxy(MR3);
kspaceR3 = getStackedKspace(MxyR3,nXCal,nYCal,nCoils);
imgR3 = getStackedImage(kSpaceR3);

R = 3;
Gamma = eye(nCoils,nCoils);
[imgSenseR3,gMapR3] = reconstructSense(sensitivityMaps,imgR3,Gamma,R);
mean(abs(gMapR3(mask)))

% visualize result
plotImage(rescale(abs(imgSenseR3.*mask)),'gray');
plotImage(medfilt2(abs(gMapR3.*mask)),'gray');
colorbar;

% R4
MR4 = load('data/r4_8_coils.mat').M;
nXCal = 256;
nYCal = 128;

MxyR4 = getMxy(MR4);
kspaceR4 = getStackedKspace(MxyR4,nXCal,nYCal,nCoils);
imgR4 = getStackedImage(kSpaceR4);

R = 4;
Gamma = eye(nCoils,nCoils);
[imgSenseR4,gMapR4] = reconstructSense(sensitivityMaps,imgR4,Gamma,R);
mean(abs(gMapR4(mask)))

% visualize result
plotImage(rescale(abs(imgSenseR4.*mask)),'gray');
plotImage(medfilt2(abs(gMapR4.*mask)),'gray');
colorbar;
