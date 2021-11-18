function [sensitivityMaps] = computeCoilSensitivities(stackImage)
    rsosImage  = rootSOSFromStacked(stackImage);
    sensitivityMaps = stackImage./rsosImage;
end

