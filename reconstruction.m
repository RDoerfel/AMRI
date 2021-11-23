%% CHANGE PATHS:
addpath(genpath('/zhome/19/2/163670/Downloads'))
addpath(genpath('/zhome/19/2/163670/myJEMEXDIR'))

%% Global Variables
clc; close all; clearvars;

nXCal = 256;
nYCal = 256;
nCoils = 4;

% Sensitivity Maps
% load('calibration_4_coils.mat');
% load('mask.mat');
load('maskCalibration.mat', 'nXCal', 'nYCal', 'mask', 'maskCalibration');

%% Calibration Acquisition:
% Mxy = getMxy(M);
% 
% [stackKspaceCal,stackImageCal] = getMatrices(Mxy,nXCal,nYCal,nCoils);
% 
% [sensitivityMaps] = computeCoilSensitivities(stackImageCal);
% imageReconstructed = lsReconstruct(stackImageCal,sensitivityMaps);
% 
% figure("Units", "Normalized","Position", [0.3, 0.2, 0.4, 0.7])
% imagesc(rescale(abs(imageReconstructed).*maskCalibration));
% colormap('gray');
% 
% figure("Units", "Normalized","Position", [0.05, 0.1, 0.4, 0.7])
% showGrid(abs(sensitivityMaps).*maskCalibration,'gray');
% sgtitle("Original Sensitivity Maps")
% 
% figure("Units", "Normalized","Position", [0.45, 0.1, 0.4, 0.7])
% sensitivityMaps = circshift(sensitivityMaps, -floor(nXCal/2), 1);
% sgtitle({' Sensitivity Maps circularly shifted along the y direction','(central row now at the top)'})
% 
% showGrid(abs(sensitivityMaps).*circshift(maskCalibration, -floor(nXCal/2), 1),'gray');
% 
% figure
% showGrid(abs(stackImageCal).*maskCalibration,'gray');
% close all


%% Noise-free vs 0.001 noise simulations (for different R's):

% Load Noise-free and noisy simulations:
% Calibration Acquisition (R=1, 256 x 256 pts)
M_noise0_R1 = load('EPI_R1_noise0.mat', 'M').M;
M_noise001_R1 = load('EPI_R1_noise001.mat', 'M').M;
% Accelerated acquision(R=2, 128 x 256 pts)
M_noise0_R2 = load('EPI_R2_noise0.mat', 'M').M; nXPI = 256; nYPI = 128;
M_noise001_R2 = load('EPI_R2_noise001.mat', 'M').M;

% Get the complex transverse magnetizations:
Mxy_noise0_R1   = getMxy(M_noise0_R1); clear M_noise0_R1
Mxy_noise001_R1 = getMxy(M_noise001_R1); clear M_noise001_R1
Mxy_noise0_R2   = getMxy(M_noise0_R2); clear M_noise0_R2
Mxy_noise001_R2 = getMxy(M_noise001_R2); clear M_noise001_R2

% Get the reduced-FOV images from all coils:
[~, CoilImages_noise0_R1]   = getMatrices(Mxy_noise0_R1,nXCal,nYCal,nCoils);
[~, CoilImages_noise001_R1] = getMatrices(Mxy_noise001_R1,nXCal,nYCal,nCoils);
[~, CoilImages_noise0_R2]   = getMatrices(Mxy_noise0_R2,nXPI,nYPI,nCoils);
[~, CoilImages_noise001_R2] = getMatrices(Mxy_noise001_R2,nXPI,nYPI,nCoils);

% % For a given R, display the coils' images with and without noise:
% figure("Units", "Normalized","Position", [0.1, 0.1, 0.4, 0.7])
% showGrid(abs(CoilImages_noise0_R2),'gray'); sgtitle("No noise")
% figure("Units", "Normalized","Position", [0.5, 0.1, 0.4, 0.7])
% showGrid(abs(CoilImages_noise001_R2),'gray'); sgtitle("0.001 Noise")

% Use the noise-free acquisition w/ R=1 and 256 x 256 pts as the 
% calibration acq. to compute the sensitivity maps:
sensitivityMaps = computeCoilSensitivities(CoilImages_noise0_R1);

% Display each coil's sensitivity map:
figure("Units", "Normalized","Position", [0.05, 0.1, 0.4, 0.7])
showGrid(abs(sensitivityMaps).*mask,'gray');
sgtitle("Sensitivity Maps")


%% Compute the noise correlation Matrix:

noise = CoilImages_noise001_R1 - CoilImages_noise0_R1;

% Display the noise-free and noise-only images from each coil:
figure("Units", "Normalized","Position", [0.1, 0.1, 0.4, 0.7])
showGrid(abs(CoilImages_noise0_R1),'gray'); sgtitle("R=1, No noise")
figure("Units", "Normalized","Position", [0.5, 0.1, 0.4, 0.7])
showGrid(abs(noise),'gray');
sgtitle("R=1, only noise")

% Initialize the noise correlation matrix:
noiseCorr = zeros(nCoils);
% For each coil pair ij:
for i = 1 : nCoils  
    for j = 1 : nCoils
        
        noise_i = noise(:,:, i); noise_i = noise_i(:);
        noise_j = noise(:,:, j); noise_j = noise_j(:);
        
        % For each pixel, multipy the noise seen by coil i with the complex
        % conjugate of the noise seen by coil j:
        noise_products = noise_i .* conj(noise_j);
        noiseCorr(i,j) = mean(noise_products(:));
        
    end
end

figure
imagesc(abs(noiseCorr)); colormap('jet'); colorbar
xticks(1:nCoils); yticks(1:nCoils); xlabel("Coil number"); ylabel("Coil number");
title("Noise Correlation Matrix")

%% Sense reconstruction for different R's:

% R = 1:
nYPI = nYCal; R = floor(nYCal/nYPI);
[imageSense_noise001_R1, g_map_noise001_R1] = ...
    sense(CoilImages_noise001_R1, sensitivityMaps, noiseCorr, R, nYCal,nXCal, nYPI);

%% R = 2:
nYPI = 128; R = floor(nYCal/nYPI);
[imageSense_noise001_R2, g_map_noise001_R2] = ...
    sense(CoilImages_noise001_R2, sensitivityMaps, noiseCorr, R, nYCal,nXCal, nYPI);

figure
subplot(1,2,1); imagesc((abs(imageSense_noise001_R2)))
colormap('gray');colorbar; axis image; title("SENSE-reconstructed image. 0.001 Noise, R=2")
subplot(1,2,2); imagesc((abs(g_map_noise001_R2)))
colormap('gray'); colorbar; axis image; title("g-factor map. 0.001 Noise, R=2")





%% Old PI acquisition:
clc
nXCal = 256; nYCal = 256;
nYPi = 128;
M = load('r2_4_coils.mat', 'M').M;
Mxy = getMxy(M);
[~, stackImagePI] = getMatrices(Mxy,nXCal,nYPi,nCoils);

% figure
% showGrid(abs(stackImagePI),'gray')
% sgtitle("Aliased images from each coil")

% SENSE Reconstruction

R = floor(nYCal/nYPI);
[imageSense, g_map] = sense(stackImagePI, sensitivityMaps, [], R, nYCal,nXCal, nYPI);

figure
imagesc((abs(imageSense)))
colormap('gray');colorbar; axis image; title("SENSE-reconstructed image. 0.001 Noise, R=2")

% hold on;
% plot(100,1,'*','MarkerSize', 10,'LineWidth',1,'color','red');
% plot(100,1+nXCal/2,'*','MarkerSize', 10,'LineWidth',1,'color','red');
% hold off;