
%%% M - which cluster does each pixel belong to
%%% vs - projection of clusters onto U . xs = vs'*(U-neuropil)
% shift back to original clusters using pixShift

% [Un Sv Vn] = svdecon(Uall'*Uall)
% VF = normc(Uall * Un) = U*S
% vs = lam * U*S
% lam * F = lam*U*S*V'*S0*U0' = lam*vs*V'*Sall*Vall

Sall = [];
Vall = zeros(4950*5,5000,'single');
iT = 0;
for nD = 1:5
    % pad V so they have equall tpts
    Sall = [Sall; dat{nD}.Sv];
    nT = size(dat{nD}.V,1);
    Vall(iT+[1:nT],[1:1000]+(nD-1)*1000) = dat{nD}.V;
    iT = iT+nT;
end
Sall = single(diag(Sall'));
Vall = Vall(1:iT,:);

lamF = vs' * Un(:,1:1000)' * Sall * Vall';


%for nD = 2:length(dat)
%    M0 = ShiftClusters(M,pixShift);
%end


