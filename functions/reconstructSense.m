function [imageSense,gMap] = reconstructSense(sensitivityMaps,stackImagePI,Gamma,R)
    [nYCal,nXCal,nCoils] = size(sensitivityMaps);
    [nYPi,nXPi,~] = size(stackImagePI);

    R = round(nYCal/nYPi); % Acceleration Factor
    shifted_CoilImages = circshift(stackImagePI, -floor(nYPi/2), 1);
    shifted_sensMaps = circshift(sensitivityMaps, -floor(nYCal/2), 1);
    
    imageSense = zeros(nYCal,nXCal);
    gMap = zeros(nYCal,nXCal);
    if nargin < 4
        Gamma = eye(nCoils,nCoils);
    end
    for y = 1:nYPi
        % loop over the entire left-right extent
        y_idx = y:nYPi:nYPi*R
        for x = 1:nXCal
            overlappingRows = zeros(R,1);
            overlappingRows(1) = y;
            for i = 2 : R
                overlappingRows(i) = overlappingRows(i-1) + floor(nYCal/R);
            end
            
            % get seinsitivity matrix
            S = shifted_sensMaps(overlappingRows,x,:);
            S = squeeze(S);
            S = transpose(S);
            % compute unfolding matrix U
            U = computeUnfolding(S,Gamma);
            % get Image values
            a = shifted_CoilImages(y,x,:);
            a = squeeze(a);
            % compute original image values
            v = U*a;
            % place into image
            imageSense([overlappingRows],x) =  v;
            % compute g values
            gMap([overlappingRows],x) =  computeG(S,Gamma); 
        end
    end
    imageSense = circshift(imageSense, +nYCal/2, 1);
    gMap = circshift(gMap, +nYCal/2, 1);
end

