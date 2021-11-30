function [sensitivityMaps] = computeCoilSensitivities(stackImage)
    rsosImage = sqrt(sum(abs(stackImage).^2,3));
    sensitivityMaps = stackImage./rsosImage;
end

