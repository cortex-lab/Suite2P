function [ops, stat, Fcell, FcellNeu]      = extractSignals(ops, m, stat)

ops.saveNeuropil = getOr(ops, 'saveNeuropil', 0);

Nk = numel(stat);

S = reshape(m.S, [], size(m.S, ndims(m.S)));

nBasis = size(S,2);

Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);

% maskNeu = ones(size(S,1), 1);

covL = [m.LtL m.LtS; m.LtS' m.StS];

% covLinv = inv(covL + 1e-4 * Ireg);
% covLinv = covLinv ./ repmat(diag(covLinv), 1, size(covLinv,2));
%% get signals  

nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F        = zeros(Nk, sum(ops.Nframes), 'single');
Fneu    = zeros(Nk, sum(ops.Nframes), 'single');
if ops.saveNeuropil
    Ntraces = zeros(nBasis, sum(ops.Nframes), 'single');
end
% S    = bsxfun(@times, S, maskNeu(:));

[Ly Lx] = size(ops.mimg1);

mimg1 = ops.mimg1(ops.yrange, ops.xrange);
% find indices of good clusters 
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
       break; 
    end
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    NT   = size(data,3);
    
    % process the data
    data = bsxfun(@minus, data, mimg1);
    data = my_conv2(data, ops.sig, [1 2]);
    data = bsxfun(@rdivide, data, m.sdmov);    
    data = single(reshape(data, [], NT));
    
    %
    Ftemp = zeros(Nk, NT, 'single');
    for k = 1:Nk
       ipix = stat(k).ipix(:)'; 
       if ~isempty(ipix)
           Ftemp(k,:) = stat(k).lam(:)' * data(ipix,:);
       end
    end
    %
    StU         = S' * data;
    Fdeconv     = covL\cat(1, Ftemp, StU);
    
    Fneu(:,ix + (1:NT))     = m.LtS * Fdeconv(1+Nk:end, :); % estimated neuropil
    F(:,ix + (1:NT))        = Fneu(:,ix + (1:NT)) + Fdeconv(1:Nk, :); % estimated ROI signal
    
    if ops.saveNeuropil
        Ntraces(:,ix + (1:NT)) = Fdeconv(1+Nk:end, :);
    end
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);
%% add the means back in to both neuropil and total
data = my_conv2(mimg1, ops.sig, [1 2]);
data = bsxfun(@rdivide, data, m.sdmov);
data = single(reshape(data, [], 1));

scalefactors = nan(numel(stat),1);
Ftemp = zeros(Nk, 1, 'single');
for k = 1:Nk
    ipix = stat(k).ipix(:)';
    if ~isempty(ipix)
        Ftemp(k,:) = stat(k).lam(:)' * data(ipix,1);
        scalefactors(k) = mean(m.sdmov(ipix));
    end
end
%
StU         = S' * data;
Fdeconv     = covL\cat(1, Ftemp, StU);

muS                     = m.LtS * Fdeconv(1+Nk:end, 1); % estimated neuropil

Fneu     = bsxfun(@plus, Fneu, muS); % estimated neuropil
F        = bsxfun(@plus, F,    muS+Fdeconv(1:Nk, 1));

Fneu     = bsxfun(@times, Fneu, scalefactors); % estimated neuropil
F        = bsxfun(@times, F,    scalefactors);

%% get activity stats
indNoNaN    = find(~ops.badframes);
ix          = cumsum(~ops.badframes) + 1;
ix          = ix(ops.badframes);
ix(ix>numel(indNoNaN)) = numel(indNoNaN);

F(:, ops.badframes)     = F(:, indNoNaN(ix));
Fneu(:, ops.badframes)  = Fneu(:, indNoNaN(ix));

dF = F - Fneu;
sd = std(dF, [], 2);

sk(:, 1) = skewness(dF, [], 2);
sk(:, 2) = sd/mean(sd); 
sk(:, 3) = (max(dF, [], 2)-median(dF, 2))./sd;
sk(:, 4) = (prctile(dF, 95, 2)-median(dF, 2))./sd;

for j = 1:numel(stat)
    stat(j).dFstat           = sk(j,:);
    stat(j).skew             = sk(j,1);
    stat(j).std              = sk(j,2);
    stat(j).maxMinusMed      = sk(j,3);
    stat(j).top5pcMinusMed   = sk(j,4);
    stat(j).blockstarts      = [0 cumsum(ops.Nframes)];
    stat(j).iplane           = ops.iplane;
end

%%
csumNframes = [0 cumsum(ops.Nframes)];
Fcell       = cell(1, length(ops.Nframes));
FcellNeu    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i}     = F(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeu{i}  = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end

if getOr(ops, 'saveNeuropil', 0)
    S = reshape(S, numel(ops.yrange), numel(ops.xrange), Nbasis);
    save(sprintf('%s/NEU_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, ops.iplane),  'ops', 'S', 'Ntraces', '-v7.3')
end
