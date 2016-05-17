function B0 = ShiftClusters(B,pixShift)

x=[1:size(pixShift,2)]; y=[1:size(pixShift,1)];
[xx0,yy0] = meshgrid(x,y);
xx = xx0 - pixShift(:,:,2);
yy = yy0 - pixShift(:,:,1);

B0 = interp2(x,y,B,xx,yy);
B0(isnan(B0(:))) = 0;

end