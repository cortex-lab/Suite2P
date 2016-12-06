function [Ff, res] = loadPROC(root, mname, datExp, npl, plane_um, pix_um, f0)

root    = fullfile(root, mname,datExp);
fdir    = dir(root);
blockstring = fdir(3).name;
root    = fullfile(root,blockstring);

clear Fc

ik = 0;
med0 = [];
igpl = [1:npl];
statall = [];
redcell = [];

for iplane = igpl
    fname = sprintf('F_%s_%s_plane%d_Nk*_proc.mat', mname, datExp, iplane);
    
    fs = dir((fullfile(root, fname)));
    fname = fs(1).name;
    load(fullfile(root, fname))
    
    ds = (ceil((npl+1)/2) - iplane)/npl;
    %
    Fcell = dat.F.Fcell;
    
    [NN NT] = size(Fcell{1});

    iscell  = logical(dat.cl.iscell);
    statall = cat(1, statall, dat.stat(iscell)');
    redcell = cat(1, redcell, dat.cl.redcell(iscell));
    
    for j = 1:length(Fcell)
        F = Fcell{j}(iscell, :);
        F = register_F(F', ds)';
        Fc{iplane,j} = F;
    end
    %
    icells = find(iscell);
    for ip = 1:length(icells)
        ik = ik+1;
        med0(ik,:) = [pix_um*dat.stat(icells(ip)).med plane_um * iplane];
    end
end
%%

Ls = cellfun(@(x) size(x,2), Fc);
Lmax = min(Ls(igpl, :), [], 1);
%
Ff  = [];
for iplane = igpl
    F0 = [];
    for j = [1:size(Fc,2)]
        F0 = cat(2, F0, Fc{iplane, j}(:, 1:Lmax(j)));
    end
    Ff = cat(1, Ff, F0);
end
Ff = Ff';

med = med0;
if npl>=9
    goodcell = remove_doubles2(Ff, med);
    %
    Ff      = Ff(:, goodcell);
    med     = med(goodcell, :);
    statall = statall(goodcell);
    redcell = redcell(goodcell);
end
%
nCum = [0 cumsum(Lmax)];
%
[NT NN] = size(Ff);

Ff = single(Ff);

%%
tic
[NT NN] = size(Ff);
parfor i = 1:NN
    Ff(:,i) = Ff(:,i) - ordfilt2(Ff(:,i), 31 * ceil(f0/3), true(201 * ceil(f0/3),1), 'symmetric');
    Ff(:,i) = Ff(:,i) - median(Ff(:,i)); 
end

res.med     = med;
res.statall = statall;
res.redcell = redcell;
%%

res.nCum = nCum;
