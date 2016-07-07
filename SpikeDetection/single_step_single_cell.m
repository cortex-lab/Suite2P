function [kernel, F1, dcell, Ffr]= single_step_single_cell(ops, Ff, F1, Fneu, Params, kernel, kerns, NT, npad, dcell)


Nbasis = size(kerns, 2);

[c, st]   = deconvL0(Params, double(F1), kernel);
fs = zeros(numel(F1), 1);
fs(st+1) = c;

X0 = conv(kernel, fs);
A = [X0(npad + [1:NT], 1) ones(NT,1) Fneu];
B = (Ff' * A) / (A'*A);
B(3) = min(B(3), ops.maxNeurop);  % force neuropil coefficient below a value

if nargout<=2
    % if only two outputs, re-estimate parameters and 
    F1([1:npad NT+npad+[1:npad]])=0;
    F1(npad + [1:NT])          = double(Ff  - A(:, [2 3])*B([2 3])');
    %         B = median(Ff(:,icell) - X0(npad + [1:NT]));    
    
    NT2 = numel(F1);
    X = zeros(NT2, Nbasis);
    for i = 1:Nbasis
        X0 = conv(kerns(:,i), fs);
        X(:,i) = X0(1:NT2);
    end
    Cv = (X' * X)/size(X,1);
    Fv = (X' * F1)/size(X,1);
    B2 = Cv \ Fv;
    kerneli = kerns * B2(1:Nbasis);
    kernel  = kerneli/sum(kerneli.^2).^.5 ;
else    
    imax        = find(c==0, 1);
    c          = c(1:imax-1);
    st         = 1+ st(1:imax-1) - npad;
    c(st<1 | st>NT)  =[];
    st(st<1 | st>NT)  =[];
    
    Ffr = zeros(size(Ff));
    dcell.c  = c;
    dcell.st = st;
    dcell.kernel = kernel;
    dcell.B   = B; % full neuropil subtraction coefficient
    Ffr(st,1)  = c; 
end
