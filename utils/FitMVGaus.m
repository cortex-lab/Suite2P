% fit 2D gaussian to cell with lam pixel weights
function params = FitMVGaus(iy,ix, lam)

% normalize pixel weigths
lam     = lam / sum(lam);

mu      = [sum(lam.*iy) sum(lam.*ix)];

xy      = bsxfun(@minus, [iy ix], mu);
xy      = bsxfun(@times, xy, sqrt(lam));

sigxy   = xy' * xy;

params.mu = mu;
params.sig = sigxy;
[evec,eval]=eig(2.5*params.sig);
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
xy = [cos(p'),sin(p')] * sqrt(eval) * evec'; % Transformation
xy = bsxfun(@plus,xy,mu);

params.xy = xy;
eval=diag(eval);

params.eval = eval;
    

