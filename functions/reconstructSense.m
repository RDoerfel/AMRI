function [imageSense,gMap] = reconstructSense(sensitivityMaps,stackImagePI,Gamma,R)
    [nYCal,nXCal,nCoils] = size(sensitivityMaps);
    imageSense = zeros(nYCal,nXCal);
    gMap = zeros(nYCal,nXCal);
    if nargin < 4
        Gamma = eye(nCoils,nCoils);
    end
    for y = 1:floor(nXCal/R)
        % loop over the entire left-right extent
        y_idx = y:floor(nXCal/R):nXCal;
        for x = 1:nXCal
            % get seinsitivity matrix
            S = sensitivityMaps([y_idx],x,:);
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
            imageSense([y_idx],x) =  v;
            % compute g values
            gMap([y_idx],x) =  computeG(S,Gamma); 
        end
    end
end

