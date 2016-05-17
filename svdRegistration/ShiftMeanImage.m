function B0 = ShiftMeanImage(A,B,pixShift,npix)


x=[1:npix]; y=[1:npix];
[xx0,yy0] = meshgrid(x,y);
xx = xx0 - pixShift(:,:,2);
yy = yy0 - pixShift(:,:,1);

B0 = interp2(x,y,B,xx,yy);
B0(isnan(B0(:))) = 0;
% subplot(1,2,1),
% imagesc(A-B)
% title('original')
% subplot(1,2,2),
% imagesc(A-B0)
% title('post-shift')
%[140 405] -> [139 406]

[cx,ix]=regZ(A,single(B));
fprintf('\ncost original: %2.4f\n',cx)
[cx,ix]=regZ(A,single(B0));
fprintf('cost post-shift: %2.4f\n',cx)

end