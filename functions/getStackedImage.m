function stackImage = getStackedImage(stackKspace)
    [nY,nX,nCoils] = size(stackKspace);
    stackImage = zeros(nY,nX,nCoils);
    for i = 1:nCoils
        stackImage(:,:,i) = kspace2image(stackKspace(:,:,i));
    end
end

