

if ~exist('loaded')
    [Uall,dat] = LoadMultiDayWhopper;
    loaded=1;
    disp('data loaded');
end

%%
addpath('~/grive/important_functions');

isGPU = 1;
nBest = 2;
npix = 512;

% shift svds to mean image of nBest and multiply Uall by Sv
CompareShiftSVDs(Uall,dat,nBest,isGPU);

U2 = Uall(ops.yrange,ops.xrange,:,:);
%% take top 1000 components of Uall
U2 = reshape(U2,[],size(U2,3)*size(U2,4));
[Un, Sv, ~] = svdecon(U2'*U2);
Sv = diag(Sv);
%%
U = U2 * Un(:, 1:1000);
clear U2;
U = normc(U);
Sv = Sv/(length(dat)/2);


%%
clf;
[ops,res] = fast_clustering_with_neuropil(ops,U,Sv);


B = single(dat{nBest}.ops.mimg1);
for nD = dinds
    A = single(dat{nD}.ops.mimg1);
    [res0,pixInv] = ShiftMasks(res,ops,A,B,isGPU);
    dat{nD}.res    = res0;
    dat{nD}.pixInv = pixInv;
end












