%% ITeRATIVE OPTIMIZATION OF PARAMETERS
addpath('D:\CODE\MariusBox\SpikeDetection')
%
% find_kernel;
%%
[NT NN] = size(Ff);
Fsp = zeros(NT, NN, 'single');
Fin = zeros(NT, NN, 'single');
Frec = zeros(NT, NN, 'single');

tic
%%
for icell = 1:20 %NN
    kernel = [0 exp(-[1:30]/6)]';
    
    kernel = kernel/sum(kernel.^2).^.5 ;
    
    F =Ff(:,icell);
    Fneu = Ffn(:,icell);
    F = F - Fneu;
    % F = F - 2.5;
    
    % F = my_conv(F', 1)';
    F = F - median(F);
    Fneu = Fneu - median(Fneu);
    sd = median(abs(F - median(F)));
    
    %
    Finp = F/sd;
    Fneu = Fneu/sd;
    Finp = double(Finp);
    kernel = double(kernel);
    
    % nt0 = numel(kernel);
    npad = 250;
    
    F1 = cat(1, zeros(npad,1), Finp, zeros(npad,1));
    Fneu = cat(1, zeros(npad,1), Fneu, zeros(npad,1));
    
    kerneli = kernel;
    F1i = F1;
    nt0 = numel(kernel);
    
    taus = [1 2 4 8 16 32];
    Nbasis = numel(taus);
    kerns = zeros(nt0, Nbasis);
    for i = 1:Nbasis
        kerns(:,i) = exp(-[1:nt0]/taus(i));
    end
    %
    NTi = numel(F1i);
    for iter = 1:3
        [F0, c, st] = deconvolveL0(F1i, kerneli, 'nonneg', 5);
        
        
        fs = zeros(numel(F1i), 1);
        fs(st+1) = c;
        
        X = zeros(NTi, Nbasis+1);
        for i = 1:Nbasis
            X0 = conv(kerns(:,i), fs);
            X(:,i) = X0(1:NTi);
        end
        
        X(:, Nbasis+1)= 1;
        X(:, Nbasis+2)= Fneu;
        
        if iter>1
            X(:, [Nbasis+1 Nbasis+2]) = [];
        end
        Cv = X' * X;
        B = Cv\(X' * F1);
        
        if iter==1
            F1i = F1 - B(Nbasis+1) - B(Nbasis+2)*Fneu;
%             B(Nbasis+2)
        end
        kerneli = kerns * B(1:Nbasis);
        kerneli = kerneli/sum(kerneli.^2).^.5 ;
    end
    
    
    F0 = F0(npad+1:end-npad);
    st = st - npad;
    c(st<1 | st>numel(F)) = [];
    st(st<1 | st>numel(F)) = [];
    
    fs = zeros(numel(F), 1);
    fs(st) = c;
    
    dcell{icell}.F0   = Finp(2:end);
    dcell{icell}.c = c;
    dcell{icell}.st = st;
    dcell{icell}.kernel = kerneli;
    
    Fsp(:, icell) = fs;
    Fr = conv(fs, kerneli);
    Frec(:, icell) = Fr(1:NT);
    Fin(:, icell) = Finp;
end
toc
%%

for icell = 1:10:6219
    plot(dcell{icell}.kernel)
    hold all
end
%%
clf

raster_transients(Fsp(900:1000, 1:1:14000)', .75)

set(gca, 'xtick', [0 50 100], 'xticklabel', [0 50 100]/2.5)

xlabel('time (seconds)')
%%



