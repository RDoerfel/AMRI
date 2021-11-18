function [image] = kspace2image(kspace)
%KSPACE2IMAGE Summary of this function goes here

    % first fftshift -> center of k-space in middle
    % second fftshift -> center (0) is in the middle
    image = fftshift(ifft2(fftshift(kspace)));

end

