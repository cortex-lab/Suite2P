%% ITeRATIVE OPTIMIZATION OF PARAMETERS
addpath('D:\CODE\MariusBox\SpikeDetection')
%
% find_kernel;
%%
tic
for icell = 1:11
    %%
%     icell = 1;
    
    kernel = kernel/sum(kernel.^2).^.5 ;
    F = dcell{icell}.F;
    % F = F - 2.5;
    
    % F = my_conv(F', 1)';
    F = F - median(F);
    sd = median(abs(F - median(F)));
    
    %
    Finp = F/sd;
    Finp = double(Finp);
    kernel = double(kernel);
    
    % nt0 = numel(kernel);
    npad = 250;
    
    F1 = cat(1, zeros(npad,1), Finp, zeros(npad,1));
    
    kerneli = kernel;
    F1i = F1;
    nt0 = numel(kernel);
    
    taus = [1 2 4 8 16 32 64];
    Nbasis = numel(taus);
    kerns = zeros(nt0, Nbasis);
    for i = 1:Nbasis
        kerns(:,i) = exp(-[1:nt0]/taus(i));
    end
    %%
    NT = numel(F1i);
    for iter = 1:3
        %%
        tic
        [F0, c, st] = deconvolveL0(F1i, kerneli, 'nonneg', 3);
        toc
        %%
        fs = zeros(numel(F1i), 1);
        fs(st+1) = c;
        
        X = zeros(NT, Nbasis+1);
        for i = 1:Nbasis
            X0 = conv(kerns(:,i), fs);
            X(:,i) = X0(1:NT);
        end        
        X(:, Nbasis+1)= 1;
        
        if iter>1
            X(:, Nbasis+1) = [];
        end
        Cv = X' * X;
        B = Cv\(X' * F1);
        
        if iter==1
            F1i = F1 - B(Nbasis+1);
        end
        kerneli = kerns * B(1:Nbasis);
        kerneli = kerneli/sum(kerneli.^2).^.5 ;
    end
    
%     figure(2)
%     plot(kerneli)
%     drawnow
%     numel(st);
    
    F0 = F0(npad+1:end-npad);
    st = st - npad;
    c(st<1 | st>numel(F)) = [];
    st(st<1 | st>numel(F)) = [];
    %
    hd = histc(dcell{icell}.tspikes, 1:1:numel(F));
    fs = zeros(numel(F), 1);
    fs(st) = c;
    
    dcell{icell}.Frec = conv(fs, kerneli);
    dcell{icell}.F0   = Finp(2:end);
    
    %
    sig = [1 3 5 7 11 15 19 25 31 41 61 121];
    cfit = zeros(length(sig), 2);
    
    for i = 1:length(sig)
        CC = corrcoef(smooth(fs(1:end-1)', sig(i) + 2), smooth(hd(2:end)', sig(i)));
        cfit(i,1) = CC(1,2).^2;
        
        CC = corrcoef(smooth(F(9:end)', sig(i) + 2), smooth(hd(1:end-8)', sig(i)));
        cfit(i,2) = CC(1,2).^2;
    end
    
    dcell{icell}.cfit = cfit;
    %
    
    s1 = fs;
    s2 = hd;
    s1 = zscore(s1, 1, 1)/numel(s1).^.5;
    s2 = zscore(s2, 1, 1)/numel(s2).^.5;
    cc = real(ifft(fft(s1, [], 1) .* fft(s2(end:-1:1), [], 1)));
    
    dcell{icell}.ccfit = cc;
    
    s1 = F;
    s2 = hd;
    s1 = zscore(s1, 1, 1)/numel(s1).^.5;
    s2 = zscore(s2, 1, 1)/numel(s2).^.5;
    cc = real(ifft(fft(s1, [], 1) .* fft(s2(end:-1:1), [], 1)));    
        
    dcell{icell}.ccfit(:,2) = cc;
    dcell{icell}.kernel = kerneli;    
    
    dcell{icell}.fs = fs;
    dcell{icell}.hd = hd;
end
toc
%%
clear h
clf

cl = colormap('jet');
cl = cl(round(linspace(1,64, numel(dcell))), :);

clf
h{5} = subplot(4,4,4+5);
for icell = 1:length(dcell)
    plot((1:1:nt0)/60, dcell{icell}.kernel)
    hold all
end
axis tight
xlabel('time (s)')
ylabel('response')

title('Kernels')
legend('cell 1', 'cell 2','cell 3', 'etc', 'Location', 'Best')
legend boxoff


h{6} = subplot(4,4,4+6);
for icell = 1:length(dcell)
    semilogx(sig/60, dcell{icell}.cfit(:,1), 'Color', cl(icell,:))
    hold on
    semilogx(sig/60, dcell{icell}.cfit(:,2), '--', 'Color', cl(icell,:))
end
axis tight

ylabel('variance explained')
xlabel('bin size (s)')
title('of ground truth')

set(gca, 'xtick', [.1 .5 1 2])

legend('model', 'raw F', 'Location', 'Best')
legend boxoff
%
h{13} = subplot(4,4,4+9);
for icell = 1:length(dcell)
    ccfit = dcell{icell}.ccfit(:,2);
    plot([-99:1:100]/60, ccfit([end-99:end 1:100]), '--', 'Color', cl(icell,:))
    hold on
end
axis tight

xlabel('time lag (s)')
ylabel('cross-correlation')

title({'cross-correlogram', 'dF'})

h{14} = subplot(4,4,4+10);

for icell = 1:length(dcell)
    ccfit = dcell{icell}.ccfit(:,1);
    plot([-99:1:100]/60, ccfit([end-99:end 1:100]), '-', 'Color', cl(icell,:))
    hold on
end
axis tight

xlabel('time lag (s)')
ylabel('cross-correlation')

title({'cross-correlogram', 'model'})


h{15} = subplot(4,4,4+11);

for icell = 1:length(dcell)
    ccfit = dcell{icell}.ccfit(:,1);
    plot((1 + [-12:1:9])/60, ccfit([end-12:end 1:9]), '-', 'Color', cl(icell,:))
    hold on
end
axis tight

xlabel('time lag (s)')
ylabel('cross-correlation')

title({'cross-correlogram', 'model'})



icell = 1;
xs = {[715*60:720*60], [1*60:80*60], [635*60:675*60], [17444:19740]};
ttls = {'bursts zoomin', 'burst zoomin', 'big dynamic range', 'failing at the limit of (structured) noise'};

for k = 1:numel(xs)
    h{k} = subplot(4,2,k);
    
    % h{1} = subplot(4,4,1);
    plot([1:length(xs{k})]/60, dcell{icell}.F0(xs{k}), 'Linewidth', .5)
    hold all
    plot([1:length(xs{k})]/60, dcell{icell}.Frec(xs{k}), 'Linewidth', 1)
    
    Fsparse = dcell{icell}.fs/mean(dcell{icell}.fs(dcell{icell}.fs>1e-4));
    plot([1:length(xs{k})]/60, 8 * Fsparse(xs{k}), 'Linewidth', 1)
    
    hd = dcell{icell}.hd(xs{k});
    for j = 1:length(hd)
        if hd(j)>0
            plot([j j]/60, 10 + [0 7]*hd(j), 'k', 'Linewidth', 1)
        end
    end
    axis tight
    if k==1
        xlabel('time (s)                                                      ')
        legend('dF', 'reconstruction', 'deconvolved', 'ground truth', 'Location', 'Best')
        legend boxoff
        
    end
    title(ttls{k})
    if k==3
       ylim([-5 20]) 
    end
end
%%


%%
ix = 0;
labl = 'ABCDEFGHIJKL';
for k = 1:length(h)
   if ~isempty(h{k})
       ix = ix+1;
       ps = h{k}.Position;
       axes('Position', [ps(1)-.05 ps(2)+ps(4) + .03, .01 .01])
       text(0, 0, labl(ix), 'Fontsize', 16)
       axis off
   end
end


