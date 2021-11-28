function fig = plotImage(data,cmap)
    fig = figure;
    imagesc(data);
    colormap(cmap);
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    axis image;
end

