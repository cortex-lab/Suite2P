function [kernel, tau, ca] = estimateKernel(ops, ca, tlag)

if nargin<2
    tlag = 1;
end

ops.fs = getOr(ops, 'fs', 30);
ops.mtau = getOr(ops, 'mtau', 1);
ops.recomputeKernel  = getOr(ops, 'recomputeKernel', 1);
ops.running_baseline = getOr(ops, 'running_baseline', 0);

if ops.running_baseline>0
    Fneu    = ca;
    NT      = size(ca,1);
    ntBase  = 2*ceil(ops.running_baseline * ops.fs/2)+1;
    Fneu    = cat(1, Fneu((ntBase-1)/2:-1:1, :), Fneu, Fneu(end:-1:end-(ntBase-1)/2, :));
    Fbase   = my_conv2(Fneu, 3, 1); % medfilt1(Fneu(:), 5);
    Fbase   = ordfilt2(Fbase, 1, true(ntBase,1));
    
    Fneu   = Fneu((ntBase-1)/2 + [1:NT], :);
    Fbase  = Fbase((ntBase-1)/2 + [1:NT], :);
    
    ca = Fneu - Fbase;
end

% ca = ca + .25 * randn(size(ca)) * std(ca(:));

if ops.recomputeKernel
    [kernel, tau] = get1Dkernel(ca, ceil(ops.mtau*ops.fs), ceil(tlag*ops.fs), ops.fs);
end




