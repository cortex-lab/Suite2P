cd('D:\DATA\F\M160825_MP027\2016-12-15\1_2_3')

st = [];
statLabels = {'iscell', 'std', 'skew', 'cmpct', 'footprint', 'mimgProjAbs'};

for iplane = 2:12    
    fname = sprintf('F_M160825_MP027_2016-12-15_plane%d.mat', iplane);
    load(fname)

    ilbl = false(1, numel(statLabels));
    st0 = [];
    for j = 1:numel(statLabels)
        if isfield(stat, statLabels{j})
            ilbl(j) = true;
            st0(:,j) = eval(sprintf('[stat.%s]', statLabels{j}));
        end
    end
    
    st = cat(1, st, st0);
end

st(isnan(st)) = 2;
% save('C:\CODE\GitHub\Suite2P\configFiles\cl_default.mat', 'st', 'statLabels')
%%

prior = buildaLDA([], st, statLabels);
prior.n = 500 * model.p;
%%
save('C:\CODE\GitHub\Suite2P\configFiles\prior_default.mat', 'prior', 'statLabels')

st = [];
save('C:\CODE\GitHub\Suite2P\configFiles\cl_default.mat', 'st', 'statLabels', 'prior')

%%
% [Ypred, ps] = evaluateLDA(model, X);