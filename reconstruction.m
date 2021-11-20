clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 128;
nCoils = 4;

%% Sensitivity Maps
load('data/calibration_4_coils.mat');
load('data/mask.mat');
load('data/maskCalibration.mat');

Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspaceCal,stackImageCal] = getMatrices(Mxy,nXCal,nYCal,nCoils);

[sensitivityMaps] = computeCoilSensitivities(stackImageCal);
imageReconstructed = lsReconstruct(stackImageCal,sensitivityMaps);

figure;
imagesc(rescale(abs(imageReconstructed).*maskCalibration));
colormap('gray');
showGrid(abs(sensitivityMaps).*maskCalibration,'jet');
showGrid(abs(stackImageCal).*maskCalibration,'gray');

%% PI
nXPi = 256;
nYPi = 128;
load('data/r2_4_coils.mat');
Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspacePI,stackImagePI] = getMatrices(Mxy,nXPi,nYPi,nCoils);

showGrid(abs(stackImagePI),'gray');

%% Sense
imageSense = zeros(nYPi,nXPi);
for y = 1:nYCal/2
    % loop over the entire left-right extent
    for x = 1:nXCal
        if(abs(squeeze(sensitivityMaps([y y+nYCal/2],x,:)))>0)
            S = transpose(squeeze(sensitivityMaps([y y+nYCal/2],x,:)));
            sol = pinv(S)*squeeze(stackImagePI(y,x,:));
            imageSense(y,x) = sol(1);
        end
    end
end
figure;
imagesc(rescale(abs(imageSense)));
