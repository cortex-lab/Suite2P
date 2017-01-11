function params = FitEllipse(ix,iy)
ix     = ix(:);
iy     = iy(:);


params = EllipseDirectFit([iy ix]);

if params(1)<0
    params = -params;
end
a = params(1);
b = params(2)/2;
c = params(3);
d = params(4)/2;
e = params(5)/2;
f = params(6);
clear params;
% checking if the ellipse is a proper real ellipse
emat = [a b d; b c e; d e f];
lmat = [a b; b c];
good_ellipse = 0;
if det(emat)~=0 && det(lmat)>0 && det(emat)/(a+c)<0
    good_ellipse = 1;
end

% find minor/major axis and center
xc = (c*d - b*e)/(b^2 - a*c);
yc = (a*e - b*d)/(b^2 - a*c);
u0 = 2*(a*e^2 + c*d^2 + f*b^2 - 2*b*d*e - a*c*f);
ra = sqrt(u0/((b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c))));
rb = sqrt(u0/((b^2-a*c)*(-1*sqrt((a-c)^2+4*b^2)-(a+c))));
% angle
if b==0 && a<c
    ang = 0;
elseif b==0 && a>c
    ang = pi/2;
elseif b~=0 && a<c
    ang = 0.5*acot((a-c)/(2*b));
else
    ang = 0.5*acot((a-c)/(2*b)) + pi/2;
end

params.ra = ra;
params.rb = rb;
params.ang = ang;
params.xc  = xc;
params.yc  = yc;
params.isgood = good_ellipse;
