function [dcell, Ffr, kernel] = run_deconvolution(Ff, f0, kernel)



Params = [1 3 3 3000]; %type, Th, Thi, maxiter

kernel = interp1(1:numel(kernel), kernel, ...
    linspace(1, numel(kernel), ceil(f0/3 * numel(kernel))));
kernel = normc(kernel');
%%
clear dcell
tic
[NT, NN] = size(Ff);
Ffr = zeros(size(Ff));
for icell = 1:NN
    F1 = Ff(:,icell);
   
    npad = 250;
    
    F1 = cat(1, zeros(npad,1), F1, zeros(npad,1));
    
    for iter = 1
        [c, st]   = deconvL0(Params, double(F1), kernel);
%         [F0, c, st] = deconvolveL0(F1i, kernel, 'nonneg', 3);
        
        fs = zeros(numel(F1), 1);
        fs(st+1) = c;
                
        X0 = conv(kernel, fs);
        B = median(Ff(:,icell) - X0(npad + [1:NT]));
    end
    
    F0          = double(Ff(:,icell)  - B);
    [c1, st1]   = deconvL0(Params, F0, kernel);
%     keyboard;
    imax    = find(c1==0, 1);
    c1      = c1(1:imax-1);
    st1     = 1+ st1(1:imax-1);
    
    dcell{icell}.c  = c1;
    dcell{icell}.st = st1;
    dcell{icell}.B  = B;
    Ffr(st1,icell)  = c1;
    
    if rem(icell, 1000)==1
        fprintf('%d neurons done \n', icell)
    end
end

dcell{1}.kernel = kernel;