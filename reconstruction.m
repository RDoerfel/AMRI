clc;
close all;
clear all;
addpath('./functions');

%% Global Variables
nXCal = 256;
nYCal = 128; % Right now: 256 (it will be overwritten by 'maskCalibration.mat')
nCoils = 4;

% Sensitivity Maps
load('data/calibration_4_coils.mat');
load('data/mask.mat');
load('data/maskCalibration.mat');

Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspaceCal,stackImageCal] = getMatrices(Mxy,nXCal,nYCal,nCoils);

[sensitivityMaps] = computeCoilSensitivities(stackImageCal);
imageReconstructed = lsReconstruct(stackImageCal,sensitivityMaps);

figure("Units", "Normalized","Position", [0.3, 0.2, 0.4, 0.7])
imagesc(rescale(abs(imageReconstructed).*maskCalibration));
colormap('gray');

figure("Units", "Normalized","Position", [0.05, 0.1, 0.4, 0.7])
showGrid(abs(sensitivityMaps).*maskCalibration,'gray');
sgtitle("Original Sensitivity Maps")

figure("Units", "Normalized","Position", [0.45, 0.1, 0.4, 0.7])
sensitivityMaps = circshift(sensitivityMaps, -floor(nXCal/2), 1);
sgtitle({' Sensitivity Maps circularly shifted along the y direction','(central row now at the top)'})

showGrid(abs(sensitivityMaps).*circshift(maskCalibration, -floor(nXCal/2), 1),'gray');
%%
figure
showGrid(abs(stackImageCal).*maskCalibration,'gray');

%% PI
nXPi = 256;
nYPi = 128;
load('data/r2_4_coils.mat');
Mx = squeeze(M(:,1,:));
My = squeeze(M(:,2,:));
Mxy = Mx + 1i*My;

[stackKspacePI,stackImagePI] = getMatrices(Mxy,nXPi,nYPi,nCoils);

figure("Units", "Normalized","Position", [0.05, 0.1, 0.4, 0.7])
showGrid(abs(stackImagePI),'gray');
sgtitle("Aliased images from each coil")

figure("Units", "Normalized","Position", [0.45, 0.1, 0.4, 0.7])
stackImagePI = circshift(stackImagePI, -floor(nYPi/2), 1);
showGrid(abs(stackImagePI),'gray');
sgtitle({'Aliased images from each coil circularly shifted along the y direction','(central row now at the top)'})

%% SENSE Reconstruction

% Initialize the full-FOV image:
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
        % get values of current pixel in each coil's aliased image:
        a = stackImagePI(y,x,:);
        a = squeeze(a);
        % compute original image values at pixels (x,y) and (x,y+FOV/2)
        v = U*a;
        % place into the full-FOV image
        imageSense([y y + nXCal/2],x) =  [v(1), v(2)];
        
    end
end

% Circularly shift back the image (central line in the image = central row in the matrix)
imageSense = circshift(imageSense, +nYCal/2, 1);
figure
imagesc(rescale(abs(imageSense))); colormap('gray'); axis image
hold on;
plot(100,1,'*','MarkerSize', 10,'LineWidth',1,'color','red');
plot(100,1+nXCal/2,'*','MarkerSize', 10,'LineWidth',1,'color','red');
hold off;