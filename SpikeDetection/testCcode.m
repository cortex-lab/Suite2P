%%
if ~exist('dcell')
    load('D:\CODE\MariusBox\BigNeuralCode\results\driv1112.mat')
    
    NT = 3e4;
    NN = numel(dcell);
    F = zeros(NT,NN, 'single');
    
    for nn = 1:NN
        F(dcell{nn}.st, nn) = dcell{nn}.c;
    end
    
    NT = find(sum(F,2)<1e-6 & [1:NT]'>100, 1) - 7;
    
    F = F(1:NT, :);
    
    % load('D:\CODE\MariusBox\BigNeuralCode\results\driv1112F.mat')
    cd D:\CODE\MariusBox\SpikeDetection
    %
    addpath('D:\CODE\MariusBox\SpikeDetection')
    Params = [1 3 3 3000]; %type, Th, Thi, maxiter
    
    rperm = randperm(NN, 20);
    for i = 1:20
        kernel = dcell{rperm(i)}.kernel;
        plot(kernel)
        hold all
    end
    % F1i(1) = 0;
    kernel = dcell{10}.kernel;
    
    [NT, NN] = size(Ff);
    dcell = cell(NN, 1);
    
    Ffr = zeros(size(Ff));
    
    tic
    for i = 1:NN
        F0          = double(Ff(:,i));
        [c1, st1]   = deconvL0(Params, F0, kernel);
        
        imax = find(c1==0, 1);
        c1 = c1(1:imax-1);
        st1 = 1+ st1(1:imax-1);
        
        dcell{i}.c  = c1;
        dcell{i}.st = st1;
        Ffr(dcell{i}.st,i) = dcell{i}.c;
        if rem(i, 1000)==1
            fprintf('%d neurons done \n', i)
        end
    end
    dcell{1}.kernel = kernel;
    toc
end
%%
NT = 3e4;
NN = numel(dcell);
F = zeros(NT,NN, 'single');

for nn = 1:NN
    F(dcell{nn}.st, nn) = dcell{nn}.c;
end

NT = find(sum(F,2)<1e-6 & [1:NT]'>100, 1);

F = F(1:NT, :);

%%
