
function redraw_figure(h)
I = hsv2rgb(cat(3, h.dat.img1.H, h.dat.img1.Sat, h.dat.img1.V));
I = min(I, 1);
axes(h.axes2); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
I = hsv2rgb(cat(3, h.dat.img2.H, h.dat.img2.Sat, h.dat.img2.V));
I = min(I, 1);
axes(h.axes3); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off