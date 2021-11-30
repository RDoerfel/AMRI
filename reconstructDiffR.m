clc;
close all;
clear all;
addpath('./functions');
addpath('./data');

%% red Sensemaps
for i = 1:8
     a = h5read ('sensmaps.h5', sprintf('/maps/magnitude/%02d',i-1));
     p = h5read ('sensmaps.h5', sprintf('/maps/phase/%02d',i-1));
     SMAP(:,:,i)=fliplr(transpose(a.*exp(sqrt(-1)*(p))));
end

%% Global Variables
nXCal = 128;
nYCal = 64;
nCoils = 8;

M = load('data/spinEcho_EPI_8ch_CAL_10T2_0CS_0T2star.mat').M;
mask = load('data/maskCalibration.mat').maskCalibration;

Mxy = getMxy(M);
stackKspace = getStackedKspace(Mxy,nXCal,nYCal,nCoils);
stackImageCal = getStackedImage(stackKspace);
stackImageCal = imresize(stackImageCal,[nXCal,2*nYCal]);

sosImage = rootSOSFromStacked(stackImageCal);
thresh = 0.05*max(abs(sosImage(:)));
mask = abs(sosImage) > thresh;

sensitivityMaps = computeCoilSensitivities(stackImageCal);

%% Different Rs 
MR2 = load('data/spinEcho_EPI_8ch_R2_10T2_0CS_0T2star.mat').M;
nXPi = 128;
nYPi = 64;
R = 2;

MxyR2 = getMxy(MR2);
kspaceR2 = getStackedKspace(MxyR2,nXPi,nYPi,nCoils);
kspaceR2 = zeroFillKspace(kspaceR2,R);

imgR2 = getStackedImage(kspaceR2);

Gamma = eye(nCoils,nCoils);
[imgSenseR2,gMapR2] = reconstructSense(sensitivityMaps,imgR2,Gamma,R);
mean(abs(gMapR2(:)))

% visualize result
plotImage(rescale(abs(imgSenseR2.*mask)),'gray');
plotImage(medfilt2(abs(gMapR2)),'hot');
h = colorbar;
ylabel(h, 'g-factor');

% R3
MR3 = load('data/spinEcho_EPI_8ch_R3_10T2_0CS_0T2star.mat').M;
I =  load('data/spinEcho_EPI_8ch_R3_10T2_0CS_0T2star.mat').I;
nXPi = I(2);
nYPi = 42;
R = 3;

MxyR3 = getMxy(MR3);
kSpaceR3 = getStackedKspace(MxyR3,nXPi,nYPi,nCoils);
kSpaceR3 = zeroFillKspace(kSpaceR3,R);

imgR3 = getStackedImage(kSpaceR3);

Gamma = eye(nCoils,nCoils);
[imgSenseR3,gMapR3] = reconstructSense(sensitivityMaps,imgR3,Gamma,R);

mean(abs(gMapR3(:)))

% visualize result
plotImage(rescale(abs(imgSenseR3.*mask)),'gray');
plotImage(medfilt2(abs(gMapR3)),'hot');
h = colorbar;
ylabel(h, 'g-factor')
% R4
MR4 = load('data/spinEcho_EPI_8ch_R4_10T2_0CS_0T2star.mat').M;
I =  load('data/spinEcho_EPI_8ch_R4_10T2_0CS_0T2star.mat').I;
nXPi= I(2);
nYPi = length(I)-1;
R = 4;

MxyR4 = getMxy(MR4);
kSpaceR4 = getStackedKspace(MxyR4,nXPi,nYPi,nCoils);
kSpaceR4 = zeroFillKspace(kSpaceR4,R);

imgR4 = getStackedImage(kSpaceR4);

Gamma = eye(nCoils,nCoils);
[imgSenseR4,gMapR4] = reconstructSense(sensitivityMaps,imgR4,Gamma,R);
mean(abs(gMapR4(:)))

% visualize result
plotImage(rescale(abs(imgSenseR4.*mask)),'gray');
plotImage(medfilt2(abs(gMapR4)),'hot');
h = colorbar;
ylabel(h, 'g-factor')
