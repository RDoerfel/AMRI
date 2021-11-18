function [kspace] = signal2kspace(data,nx,ny)
%SIGNAL2MATRIX reshape from signal to kspace matrix
    kspace = reshape(data,nx,ny);
    kspace = transpose(kspace);
    kspace(2:2:end,:) = fliplr(kspace(2:2:end,:));
end

