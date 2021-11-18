clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nX = 256;
nY = 128;
nCoils = 4;

%% Sensitivity Maps
load('calibration_4_coils.mat');
load('mask.mat');

Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspaceCal,stackImageCal] = getMatrices(Mxy,nX,nY,nCoils);

[sensitivityMaps] = computeCoilSensitivities(stackImageCal);
imageReconstructed = lsReconstruct(stackImageCal,sensitivityMaps);

figure;
imagesc(abs(imageReconstructed).*mask);
showGrid(abs(sensitivityMaps).*mask);
showGrid(abs(stackImageCal).*mask);

%% reconstruct image


%% PI
load('r2_4_coils.mat');
Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspacePI,stackImagePI] = getMatrices(Mxy,nX,nY,nCoils);

showGrid(abs(stackImagePI));

%% Sense
imageSense = zeros(nY,nX);

% for x = 1:nX/2
%     for y = 1:nY
for x = 1:nX/2
    % loop over the entire left-right extent
    for y = 1:nY
        S = transpose(reshape(sensitivityMaps(y,[x x+nX/2],:),2,[]));
        imageSense(y,[x x+nX/2]) = pinv(S)*reshape(stackImagePI(z,x,1,:),[],1);
    end
end

figure;
imagesc(imageSense);
