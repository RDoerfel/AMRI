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

showGrid(abs(sensitivityMaps).*maskCalibration,'gray');
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
imageSense = zeros(nYCal,nXCal);
for y = 1:nYPi
    % loop over the entire left-right extent
    for x = 1:nXCal
        % get seinsitivity matrix
        S = sensitivityMaps([y y + nXCal/2],x,:);
        S = squeeze(S);
        S = transpose(S);
        % compute unfolding matrix U
        U = pinv(S);
        % get Image values
        a = stackImagePI(y,x,:);
        a = squeeze(a);
        % compute original image values
        v = U*a;
        % place into image
        imageSense([y y + nXCal/2],x) =  [v(1),v(2)];
    end
end
figure;
imagesc(rescale(abs(imageSense)));
hold on;
plot(100,1,'*','MarkerSize', 10,'LineWidth',1,'color','red');
plot(100,1+nXCal/2,'*','MarkerSize', 10,'LineWidth',1,'color','red');
hold off;