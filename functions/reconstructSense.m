function [imageSense,gMap] = reconstructSense(sensitivityMaps,stackImagePI,Gamma)
    [nYCal,nXCal,nCoils] = size(sensitivityMaps);
    imageSense = zeros(nYCal,nXCal);
    gMap = zeros(nYCal,nXCal);
    if nargin < 3
        Gamma = eye(nCoils,nCoils);
    end
    for y = 1:nXCal/2
        % loop over the entire left-right extent
        for x = 1:nXCal
            % get seinsitivity matrix
            S = sensitivityMaps([y y + nXCal/2],x,:);
            S = squeeze(S);
            S = transpose(S);
            % compute unfolding matrix U
            U = computeUnfolding(S,Gamma);
            % get Image values
            a = stackImagePI(y,x,:);
            a = squeeze(a);
            % compute original image values
            v = U*a;
            % place into image
            imageSense([y y + nXCal/2],x) =  v;
            % compute g values
            gMap([y y + nXCal/2],x) =  computeG(S,Gamma); 
        end
    end
end

