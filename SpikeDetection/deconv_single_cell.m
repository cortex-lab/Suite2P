function [dcell, Ffr]= deconv_single_cell(Ff, Fneu, Params, kernel)

npad = 250;
NT = size(Ff,1);
F1 = cat(1, zeros(npad,1), double(Ff), zeros(npad,1));
F1 = F1 - median(F1);

nt0 = numel(kernel);
taus = [1 2 4 8 16 32 64];
Nbasis = numel(taus);
kerns = zeros(nt0, Nbasis);
for i = 1:Nbasis
    kerns(:,i) = exp(-[1:nt0]/taus(i));
end

for iter = 1:10
    [c, st]   = deconvL0(Params, double(F1), kernel);
    %         [F0, c, st] = deconvolveL0(F1i, kernel, 'nonneg', 3);
    
    fs = zeros(numel(F1), 1);
    fs(st+1) = c;
        
    X0 = conv(kernel, fs);
    
    A = [X0(npad + [1:NT], 1) ones(NT,1) Fneu];
    B = (Ff' * A) / (A'*A);
    
    F1([1:npad NT+npad+[1:npad]])=0;
    F1(npad + [1:NT])          = double(Ff  - A(:, [2 3])*B([2 3])');
    %         B = median(Ff(:,icell) - X0(npad + [1:NT]));    
    
    NT2 = numel(F1);
    X = zeros(NT2, Nbasis);
    for i = 1:Nbasis
        X0 = conv(kerns(:,i), fs);
        X(:,i) = X0(1:NT2);
    end
    Cv = X' * X;
    B2 = Cv\(X' * F1);
    kerneli = kerns * B2(1:Nbasis);
    kernel  = kerneli/sum(kerneli.^2).^.5 ;
end

err = F1(npad + [1:NT]) - A * B';
errk = mean(err(:).^2).^.5;


[c1, st1]   = deconvL0(Params, F1, kernel);
imax        = find(c1==0, 1);
c1          = c1(1:imax-1);
st1         = 1+ st1(1:imax-1) - npad;
c1(st1<1 | st1>NT)  =[];
st1(st1<1 | st1>NT)  =[];


Ffr = zeros(size(Ff));
dcell.c  = c1;
dcell.st = st1;
dcell.kernel = kernel;
dcell.B  = B;
Ffr(st1,1)  = c1;