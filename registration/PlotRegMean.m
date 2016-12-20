function PlotRegMean(ops1,ops)

figure('position', [900 50 900 900])
    ax = ceil(sqrt(numel(ops1)));
    for i = 1:length(ops1)
        subplot(ax,ax,i)
        imagesc(ops1{i}.mimg)
        colormap('gray')
        title(sprintf('Registration for plane %d, mouse %s, date %s', ...
            i, ops.mouse_name, ops.date))
    end
    drawnow