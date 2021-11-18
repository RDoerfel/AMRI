function showGrid(gridData)

    nCoils = size(gridData,3);
    figure;
    nRows = ceil(sqrt(nCoils));
    nCols = ceil(sqrt(nCoils));
    for i = 1:nCoils
        subplot(nRows,nCols,i);
        imagesc(gridData(:,:,i));
        colormap jet;
    end
    
end

