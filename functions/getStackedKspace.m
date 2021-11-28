function [stackKspace] = getStackedKspace(Mxy,nX,nY,nCoils)
    stackKspace = zeros(nY,nX,nCoils);
    for i = 1:nCoils
        stackKspace(:,:,i) = signal2kspace(Mxy(:,i),nX,nY);
    end
end

