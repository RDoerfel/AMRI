function [rsosImage] = rootSOSFromStacked(stackedData)
    rsosImage  = sqrt(sum(abs(stackedData).^2,3));
end

