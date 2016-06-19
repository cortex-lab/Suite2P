
function redraw_meanimg(h)

if h.dat.procmap
    I = h.dat.mimg_proc(:,:,h.dat.map);
else    
    I = h.dat.mimg(:,:,h.dat.map);
end

mu = median(I(:));
sd1 = mean(abs(I(I<mu) - mu));
sd2 = mean(abs(I(I>mu) - mu));

axes(h.axes2); imagesc(I, mu + 5*[-sd1 sd2]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')
axes(h.axes3); imagesc(I, mu + 5*[-sd1 sd2]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')

drawnow