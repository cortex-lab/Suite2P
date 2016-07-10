function [FfA, AtA, F1, dcell]= ...
    single_step_single_cell(F1, Params, kernel, kerns,NT, npad, dcell)

znpad = zeros(npad,1);
F1copy = cat(1, znpad, F1, znpad); 
    
Nbasis = size(kerns,2);
F1copy = double(F1copy);

[c, st]     = deconvL0(Params, F1copy, kernel);
imax        = find(c==0, 1);
c           = c(1:imax-1);
st          = st(1:imax-1);
        
FfA = [];
AtA = [];

if nargout<=3
    % if only two outputs, re-estimate parameters
    
    fs = zeros(numel(F1copy), 1);
    fs(st+1) = c;
    NT2 = numel(F1copy);
    X = zeros(NT2, Nbasis);
    for i = 1:Nbasis
        X0 = conv(kerns(:,i), fs);
        X(:,i) = X0(1:NT2);
    end
    
    % X0 = conv(kernel, fs);
    A = X(npad + [1:NT], :) ;
    FfA = (F1' * A)/size(F1,1);
    AtA = (A'*A)/size(F1,1) + 1e-6;
    F1 = fs(1:numel(F1));
    
    %     B =  FfA/AtA ;
%     B(2) = min(B(2), maxNeurop);  % force neuropil coefficient below a value
%     B2 = B(3:end)';
%     Bout = B(1:2);
%    
%     kerneli = kerns * B2(1:Nbasis);
%     kernel = kerneli;
%     kernel  = normc(kernel) ; % renormalize the kernel
else   
    st         = st+1 - npad;
    c(st<1 | st>NT)  =[];
    st(st<1 | st>NT)  =[];
    dcell.c  = c;
    dcell.st = st;
    dcell.kernel = kernel;
end

% err = std((Ff - A * B'));
1;
