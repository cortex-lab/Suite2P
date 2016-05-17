function [cx,ixS] = SubPixel2D(cc,ixAll,nWind)
    
[nX nY] = size(cc);
nA      = [nX nY];

ixAll      = InvRegInds(ixAll,nX,nY);
ixAllS     = FftShift2D(ixAll,nX,nY,1);

% take only window around mean shift
cc0        = fftshift(cc);
dWind      = [-nWind:nWind];
xWind      = unique(min(nX,max(1,ixAllS(1)+dWind)));
yWind      = unique(min(nX,max(1,ixAllS(2)+dWind)));
cc0        = cc0(xWind,yWind);
[nXw nYw]  = size(cc0);

[cmax0,ixmax0] = max(cc0(:));
[ix1,ix2] = ind2sub(size(cc0),ixmax0);
ix0 = [ix1 ix2];

dX = [-1:1]; dY = [-1:1];
xinds = (unique(min(nXw,max(1,ix0(1)+dX))));
yinds = (unique(min(nYw,max(1,ix0(2)+dY))));

cZ=cc0(xinds,yinds);

xrange = xinds-ix0(1); 
yrange = yinds-ix0(2);
[dX2,dY2] = ndgrid(linspace(min(xrange),max(xrange),30), ...
                   linspace(min(yrange),max(yrange),30));
[dX,dY]   = ndgrid(xinds-ix0(1),yinds-ix0(2));

cZi=interp2(dY,dX,cZ,dY2,dX2,'spline');

[cx,imax] = max(cZi(:));
X2 = dX2(imax);
Y2 = dY2(imax);

% convert ix0 to real coordinates
ix0(1)    = ix0(1) + min(xWind) - 1;
ix0(2)    = ix0(2) + min(yWind) - 1;
ixS       = FftShift2D(ix0,nX,nY,-1);

ixS(1) = ixS(1) + X2;
ixS(2) = ixS(2) + Y2;

