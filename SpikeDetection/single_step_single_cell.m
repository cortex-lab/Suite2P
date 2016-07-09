function [kernel, B, err, dcell]= ...
    single_step_single_cell(maxNeurop, Ff, B, Fneu, Params, kernel, kerns, NT, npad, dcell)

F1 = zeros(NT+2*npad, 1, 'single');
F1([1:npad NT+npad+[1:npad]], 1)=0;
A(:, [2 3]) = [ones(NT,1) Fneu];
F1(npad + [1:NT])          = double(Ff  - A(:, [2 3])*B([2 3]));
    
F1copy = double(F1);

Nbasis = size(kerns, 2);

[c, st]     = deconvL0(Params, F1copy, kernel);
imax        = find(c==0, 1);
c           = c(1:imax-1);
st = st(1:imax-1);

fs = zeros(numel(F1), 1);
fs(st+1) = c;

X0 = conv(kernel, fs);
A(:,1) = X0(npad + [1:NT], 1) ;
FfA = (Ff' * A)/size(Ff,1);
AtA = (A'*A + 1e-6)/size(Ff,1);
B =  FfA/AtA ;
B(3) = min(B(3), maxNeurop);  % force neuropil coefficient below a value

if nargout<=3
    % if only two outputs, re-estimate parameters
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
    kernel = kerneli;
%     kernel = max(0, kernel); % make sure kernel is positive
            
    kernel  = normc(kernel) ; % renormalize the kernel
else   
    st         = 1+ st - npad;
    c(st<1 | st>NT)  =[];
    st(st<1 | st>NT)  =[];
    dcell.c  = c;
    dcell.st = st;
    dcell.kernel = kernel;
    dcell.B   = B; % full neuropil subtraction coefficient
end

err = sum((Ff - A * B').^2);
1;
