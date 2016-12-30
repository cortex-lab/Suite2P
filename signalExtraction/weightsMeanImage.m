function [stat, we] = weightsMeanImage(ops, stat, m)

S = reshape(m.S, [], size(m.S, ndims(m.S)));

nBasis = size(S,2);
Nk = numel(stat);
Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);

covL = [m.LtL m.LtS; m.LtS' m.StS] + 1e-4 * Ireg;
%%
mimg1 = ops.mimg1(ops.yrange, ops.xrange);
mimg1 = my_conv2(mimg1, ops.sig, [1 2]);
mimg1 = mimg1 - my_conv2(mimg1, 5*ops.diameter, [1 2]);
mimg1 = mimg1 ./ my_conv2(mimg1.^2, 5*ops.diameter, [1 2]).^.5;
% mimg1 = bsxfun(@rdivide, mimg1, m.sdmov);

mimg1 = mimg1(:);
Ftemp = zeros(Nk, 1, 'single');
for k = 1:Nk
    ipix = stat(k).ipix(:)';
    if ~isempty(ipix)
        Ftemp(k) = stat(k).lam(:)' * mimg1(ipix);
    end
end
%
StU         = S' * mimg1(:);
we          = covL \ cat(1, Ftemp, StU);

we = we(1:Nk);

for j = 1:Nk
    stat(j).mimgProj = we(j);
    % take absolute value: some cells mean activity is below neuropil!
    stat(j).mimgProjAbs = abs(we(j));
end
1;