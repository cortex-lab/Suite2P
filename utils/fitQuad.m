% fits quadratic to y = f(x) and outputs yq = f(xq)
% x = nsamps x 1
% y = nsamps x nt
% xq = nsampsq x 1
% yq = nsampsq x nt
function yq = fitQuad(x,y,xq)

X    = [ones(numel(x),1) x x.^2];
XT   = (X'*X)\X';

a    = XT * y;

XQ   = [ones(numel(xq),1) xq xq.^2];
yq   = XQ * a;