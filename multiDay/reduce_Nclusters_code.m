xs(ifeat + nFeat * [0:1:Npix-1]) = -Inf;
[M2, ifeat] = max(xs,[],1);

dM = M - M2;


icl         = get_cl(iclustNEW, Nk);

dMx = cellfun(@(x) sum(dM(x)), icl);
ix = abs(dMx - min(dMx(iSortFeat), [], 1)') < 1e-10;
ilow = true(Nk, 1);
ilow(iSortFeat(2:nFeat, ix)) = false;
ix = ix & ilow;

indrem = Nkiter(k) - Nkiter(k+1);
[dSort, isort] = sort(dMx./ix, 'descend');

indrem = min(sum(ix), indrem);
% replace each pixel in destroyed clusters with second best cluster
for j = (Nk-indrem+1):Nk
    iC = isort(j);
    iBestSecond = ifeat(icl{iC}) + nFeat * (iclust(icl{iC})-1);
    
    igood = ~ismember( iSortFeat(iBestSecond), isort((Nk-indrem+1):Nk));
    
    iclustNEW(icl{iC}(igood)) = iSortFeat(iBestSecond(igood));
    iclustNEW(icl{iC}(~igood)) = iclust(icl{iC}(~igood));
    
    if sum(ismember( iclustNEW(icl{iC}), isort((Nk-indrem+1):Nk)))>0
        keyboard;
    end
end


% compute backwards transformation
iback = zeros(1, Nk);
iback(isort) = 1:Nk;

iclust = iback(iclustNEW);

Nk = Nk - indrem;
icl         = get_cl(iclust, Nk);
iSortFeat   = get_isortfeat(icl, Nk, Ly, nFeat);