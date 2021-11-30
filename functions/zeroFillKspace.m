function [stackKspacePI] = zeroFillKspace(stackKspaceUndersampled,R)
    [nY,nX,nCoils] = size(stackKspaceUndersampled);
    nY = R*nY;
    stackKspacePI = zeros(nY,nX,nCoils);
    stackKspacePI(1:R:end,:,:) = stackKspaceUndersampled;
end

