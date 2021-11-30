function [smoothed] = smootheSensitityMaps(maps,kernelSize)
    kernel = ones(kernelSize)/kernelSize^2;
    [nx,ny,nCoils] = size(maps);
    smoothed = zeros(nx,ny,nCoils);
    for i = 1:nCoils
        smoothed(:,:,i) = conv2(maps(:,:,i),kernel,'same'); 
    end
end

