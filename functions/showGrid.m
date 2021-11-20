function showGrid(gridData,cmap)

    nCoils = size(gridData,3);
    figure;
    nRows = ceil(sqrt(nCoils));
    nCols = ceil(sqrt(nCoils));
    for i = 1:nCoils
        [nx,ny] = size(gridData(:,:,i));
        subplot(nRows,nCols,i);
        imagesc(rescale(gridData(:,:,i)));
        colormap(cmap);
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        axis image;
        
        hold on;
        plot(100,1,'*','MarkerSize', 10,'LineWidth',1,'color','red');
        if(nx > 128)
            plot(100,1+nx/2,'*','MarkerSize', 10,'LineWidth',1,'color','red');
        end
        hold off;

    end
    
end

