function [snrMap] = computeSNRMap(gMap,R)
    snrMap = SNR./(gMap.*sqrt(R));
end

