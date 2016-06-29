function redraw_fluorescence(h)
axes(h.axes4)
hold off
[NN NT] = size(h.dat.F.trace);
plot(my_conv_local(medfilt1(double(h.dat.F.trace(h.dat.F.ichosen,:)), 3), 3))
axis tight
hold on

if h.dat.plot_neu
    if isfield(h.dat.F, 'neurop')
        plot(my_conv_local(medfilt1(double(h.dat.F.neurop(h.dat.F.ichosen,:)), 3), 3))
    end
end

% plot([0 NT], [0 0], 'k', 'Linewidth', 2)
% axis off
