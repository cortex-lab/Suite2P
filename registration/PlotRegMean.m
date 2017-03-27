function PlotRegMean(ops1,ops)

figure('position', [900 50 900 900])
ax = ceil(sqrt(numel(ops1)));
for i = 1:length(ops1)
    subplot(ax,ax,i)
    imagesc(ops1{i}.mimg)
    colormap('gray')
    title(sprintf('Plane %d', i))
end
annotation('textbox', [0 .96 1 .03], 'String', ...
    sprintf('Registration for mouse %s, date %s', ops.mouse_name, ops.date), ...
    'Interpreter', 'none', 'FontWeight', 'bold', 'LineStyle', 'none', ...
    'HorizontalAlignment', 'center');
drawnow