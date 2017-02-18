function redraw_fluorescence(h)
axes(h.axes4)
hold off

ichosen = h.dat.F.ichosen;
F = [];
Fneu = [];
for j = 1:numel(h.dat.Fcell)
    F    = cat(2, F, h.dat.Fcell{j}(ichosen, :));
    Fneu = cat(2, Fneu, h.dat.FcellNeu{j}(ichosen, :));
end

plot(my_conv_local(medfilt1(double(F), 3), 3))
axis tight
hold on

coefNeu = 0.7;
if isfield(h.dat.stat, 'neuropilCoefficient')
    coefNeu = h.dat.stat(ichosen).neuropilCoefficient;
end

if isfield(h.dat, 'FcellNeu')
    plot(coefNeu * my_conv_local(medfilt1(double(Fneu), 3), 3))
end

box off
set(gca, 'xcolor', 'w')
% plot([0 NT], [0 0], 'k', 'Linewidth', 2)
% axis off
