function [imageSense, g_map] = sense(CoilImages, sensitivityMaps, noiseCorr, nYCal, nXCal, nYPI, nXPI)

    R = floor(nYCal/nYPI); % Acceleration Factor
    shifted_CoilImages = circshift(CoilImages, -floor(nYPI/2), 1);
    shifted_sensMaps = circshift(sensitivityMaps, -floor(nYCal/2), 1);
    nCoils = size(sensitivityMaps,3);

    imageSense = zeros(nYCal,nXCal); % Initialize the full-FOV image

    if ~isempty(noiseCorr) 
        % If 'noiseCorr' is provided, g-maps are to be computed:
        compute_g = 1;
    else 
        % Otherwise, assume 'noiseCorr' equals the identity matrix:  
        noiseCorr = eye(nCoils);
        compute_g = 0;
    end
    g_map = zeros(nYCal,nXCal);

    % Go through all rows of the coils' reduced-FOV images:
    for y = 1:nYPI
        % loop over the entire left-right extent
        for x = 1:nXCal
            % Identify the overlapping rows in the reduced-FOV images:
            overlappingRows = zeros(R,1);
            overlappingRows(1) = y;
            for i = 2 : R
                overlappingRows(i) = overlappingRows(i-1) + floor(nYCal/R);
            end

            % Get coils' sensitivity values at the aliased pixels:
            % S = shifted_sensitivityMaps([y y + nXCal/2],x,:);
            S = shifted_sensMaps(overlappingRows,x,:);
            S = squeeze(S);
            if size(S,1) ~= nCoils
                S = transpose(S);
            end
            S_H = ctranspose(S);

            
            % Compute g-value at the aliased pixels
            SH_invPsi_S = S_H * inv(noiseCorr) * S;
            inv_SH_invPsi_S = inv(SH_invPsi_S);

            %g_map(y,x) = sqrt( inv_SH_invPsi_S(1,1)*SH_invPsi_S(1,1) );
            %g_map(y+nXCal/2,x) = sqrt(  inv_SH_invPsi_S(2,2)*SH_invPsi_S(2,2) );
            g_map(overlappingRows,x) = sqrt( diag(inv_SH_invPsi_S).*diag(SH_invPsi_S) );


            % Compute unfolding matrix U
            % U = pinv(S);
            U = ( S_H * inv(noiseCorr) * S ) \ ( S_H * inv(noiseCorr) );

            % get values of current pixel in each coil's aliased image:
            a = shifted_CoilImages(y,x,:);
            a = squeeze(a);
            % compute original image values at the pixels that became overlapped:
            v = U*a;

            % place into the full-FOV image
            imageSense(overlappingRows,x) =  v;
        end
    end

    % Circularly shift back the imageS (central line in the image = central row in the matrix)
    imageSense = circshift(imageSense, +nYCal/2, 1);
    g_map = circshift(g_map, +nYCal/2, 1);


end

