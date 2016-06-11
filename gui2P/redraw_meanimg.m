
function redraw_meanimg(h)

I = h.dat.mimg(:,:,h.dat.map);

axes(h.axes2); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')
axes(h.axes3); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')