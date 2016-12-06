function [ fs, dcell]= ...
    single_step_single_cell(F1, Params, kernel,NT, npad, dcell)

znpad = zeros(npad,1);
F1copy = cat(1, znpad, F1, znpad); 
    
F1copy = double(F1copy);

[c, st]     = deconvL0(Params, F1copy, kernel);
imax        = find(c==0, 1);
c           = c(1:imax-1);
st          = st(1:imax-1);

st = st + 1- npad;

c(st<1 | st>NT) = [];
st(st<1 | st>NT) = [];

fs = zeros(NT, 1);
fs(st) = c;

if nargout>1
     
    dcell.c  = c;
    dcell.st = st;
    dcell.kernel = kernel;
end

% err = std((Ff - X * B'));
1;
