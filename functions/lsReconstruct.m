function [reconstructedImage] = lsReconstruct(stackImage,sensitivityMaps)
    % c = Sm -> m = S⁻¹c = (S* .* S)⁻¹.* S*
    reconstructedImage = sum(stackImage.*conj(sensitivityMaps),4)./sum(sensitivityMaps.*conj(sensitivityMaps),4);
    reconstructedImage = sum(reconstructedImage,3);
end

