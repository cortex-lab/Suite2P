function [FfA, AtA, fs, dcell]= ...
    single_step_single_cell(F1, Params, kernel, kerns,NT, npad, dcell)

znpad = zeros(npad,1);
F1copy = cat(1, znpad, F1, znpad); 
    
Nbasis = size(kerns,2);
F1copy = double(F1copy);

[c, st]     = deconvL0(Params, F1copy, kernel);
imax        = find(c==0, 1);
c           = c(1:imax-1);
st          = st(1:imax-1);

st = st + 1- npad;

c(st<1 | st>NT) = [];
st(st<1 | st>NT) = [];

FfA = [];
AtA = [];
fs = [];

if nargout<=3
    % if only two outputs, re-estimate parameters
    
    fs = zeros(NT, 1);
    fs(st) = c;
    X = zeros(NT, Nbasis);
    for i = 1:Nbasis
        X0 = conv(kerns(:,i), fs);
        X(:,i) = X0(1:NT);
    end
    
    % X0 = conv(kernel, fs);
    
    FfA = (F1' * X)/size(F1,1);
    AtA = (X'*X)/size(F1,1) + 1e-6;
    
    %     B =  FfA/AtA ;
%     B(2) = min(B(2), maxNeurop);  % force neuropil coefficient below a value
%     B2 = B(3:end)';
%     Bout = B(1:2);
%    
%     kerneli = kerns * B2(1:Nbasis);
%     kernel = kerneli;
%     kernel  = normc(kernel) ; % renormalize the kernel
else   
    dcell.c  = c;
    dcell.st = st;
    dcell.kernel = kernel;
end

% err = std((Ff - X * B'));
1;
