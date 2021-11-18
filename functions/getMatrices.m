function [stackKspace,stackImage] = getMatrices(Mxy,nX,nY,nCoils)
    stackKspace = zeros(nY,nX,nCoils);
    stackImage = zeros(nY,nX,nCoils);
    for i = 1:nCoils
        stackKspace(:,:,i) = signal2kspace(Mxy(:,i),nX,nY);
        stackImage(:,:,i) = kspace2image(stackKspace(:,:,i));
    end
end

