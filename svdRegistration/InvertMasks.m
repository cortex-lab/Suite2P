

yrange2=[10:492];
xrange2=[24:505];
nX = length(yrange2);
nY = length(xrange2);

for nD = 2%dinds
    pixShift = dat{nD}.pixShift;
    %x = [1:npix^2]'; y = [1:npix^2]';

    [x,y] = ndgrid(yrange2,xrange2);
    x=x(:); y=y(:);

    %pS = zeros(2*npix,2*npix,'single');
    %pS(1:npix,1:npix) = pixShift(:,:,1);
    %pS(npix+[1:npix],npix+[1:npix]) = pixShift(:,:,2);
    pS = [reshape(pixShift(yrange2,xrange2,1),[],1); reshape(pixShift(yrange2,xrange2,2),[],1)];

    xy = [x;y];
    xy0 = pS + xy;
    x0 = xy0(1:length(x));     y0 = xy0(length(x)+[1:length(x)]);

    pixInv = zeros(nX,nY,2,'single');
    pixInv(:,:,1) = reshape(x - x0,nX,nY);
    pixInv(:,:,2) = reshape(y - y0,nX,nY);



    %[xx,yy] = ndgrid([1:npix],[1:npix]);
    %x0 = pixShift(:,:,1).*xx;
    %y0 = pixShift(:,:,2).*yy;
    
end