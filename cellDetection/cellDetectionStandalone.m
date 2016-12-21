function [stat, F, Fneu] = cellDetectionStandalone(dat, ops)
[Ly, Lx, NT] = size(dat);

% the following options are described in the main pipeline
ops.clustrules.diameter     = getOr(ops, 'diameter', 10); % avg cell diameter
ops.sig                     = getOr(ops, 'sig', 0.5); % smoothing constant in pixels. 
ops.ResultsSavePath         = getOr(ops, 'ResultsSavePath', []); % where to save the results
ops.useGPU                  = getOr(ops, 'useGPU', 0); % whether to use a GPU
ops.fs                      = getOr(ops, 'fs', 10); % sampling rate of the movie
ops.NavgFramesSVD           = getOr(ops, 'NavgFramesSVD', max(1000, ceil(NT/ops.fs))); % size of binned movie
ops.NavgFramesSVD           = min(NT, ops.NavgFramesSVD); 
ops.fig                     = getOr(ops, 'fig', 0); % whether to show figures

navg = round(NT/ops.NavgFramesSVD);
batch = navg * ceil(2000/navg);

mov = zeros(Ly, Lx, ceil(NT/navg), 'single');
ix = 0;
iy = 0;
while ix+navg<=NT
    irange = ix + [1:batch];
    irange(irange>NT) = [];
    irange = irange(1:floor(numel(irange)/navg)*navg);
    
    mov(:,:,iy + [1:numel(irange)/navg]) = ...
        sq(mean(reshape(dat(:,:,irange), Ly, Lx, navg, []),3));
    
    ix = ix + numel(irange);
    iy = iy + numel(irange)/navg;
end
mov(:,:,iy+1:end) = [];

if ops.sig>0.05
    for i = 1:size(mov,3)
        I = mov(:,:,i);
        I = my_conv2(I, ops.sig, [1 2]);
        mov(:,:,i) = I;
    end
end

% filter temporally
mov = reshape(mov, [], size(mov,3));
nBatch = 2^15;
for j = 1:ceil(size(mov,1)/nBatch)
    ind = ((j-1)*nBatch+1):(j*nBatch);
    ind(ind>size(mov,1)) = [];
    mov(ind,:) = mov(ind,:) - my_conv2(mov(ind,:), (3*ops.fs), 2);
end

sdmov           = mean(mov.^2,2).^.5;
mov             = bsxfun(@rdivide, mov, sdmov);

model.sdmov     = reshape(sdmov, numel(ops.yrange), numel(ops.xrange));
if ops.useGPU
    COV             = gpuBlockXtX(mov)/size(mov,1);
else
    COV             = mov' * mov/size(mov,1);
end

ops.nSVDforROI = min(size(COV,1)-2, ops.nSVDforROI);

if ops.useGPU && size(COV,1)<1.2e4
    [V, ~, ~]      = svd(gpuArray(double(COV)));
    V               = single(V(:, 1:ops.nSVDforROI));
    V = gather(V);
else
    [V, ~]          = eigs(double(COV), ops.nSVDforROI);
end

if ops.useGPU
    U               = gpuBlockXY(mov, V);
else
    U               = mov * V;
end
U               = single(U);

[ops, stat, model] = sourcery(ops, U, model);

ops.clustrules      = get_clustrules(ops.clustrules);
stat                = classifyROI(stat, ops.clustrules);

% extract signals
[ops, stat,F,Fneu] = extractSigStandalone(dat, ops, model, stat);


if ~isempty(ops.ResultsSavePath)
   % save results file for GUI
   save(sprintf('%s_results.mat', ops.ResultsSavePath),  'ops',  'stat',...
            'F', 'Fneu', '-v7.3')   
end
end

function [ops, stat,F,Fneu] = extractSigStandalone(I, ops, m, stat)
% given a clustering model, extract the signals

Nk = numel(stat);
[Ly, Lx, NT] = size(I);

S = reshape(m.S, [], size(m.S, ndims(m.S)));

nBasis = size(S,2);

Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);

covL = [m.LtL m.LtS; m.LtS' m.StS];

covLinv = inv(covL + 1e-4 * Ireg);
covLinv = covLinv ./ repmat(diag(covLinv), 1, size(covLinv,2));

% get signals  
nimgbatch = 2000;

ix = 0;

tic
F        = zeros(Nk, NT, 'single');
Fneu    = zeros(Nk, NT, 'single');

while 1
    irange = ix +[1:nimgbatch];
    irange(irange>NT) = [];
    if isempty(irange)
        break;
    end
    
    data = I(:, :, irange);    
    data = single(data);
    NTbatch = size(data,3);
    
    % process the data
    data = my_conv2(data, ops.sig, [1 2]);
    data = bsxfun(@rdivide, data, m.sdmov);    
    data = single(reshape(data, [], NTbatch));
    
    %
    Ftemp = zeros(Nk, NTbatch, 'single');
    for k = 1:Nk
       ipix = stat(k).ipix(:)'; 
       if ~isempty(ipix)
           Ftemp(k,:) = stat(k).lam(:)' * data(ipix,:);
       end
    end
    %
    StU         = S' * data;
    Fdeconv     = covLinv * cat(1, Ftemp, StU);
    
    Fneu(:,irange)     = m.LtS * Fdeconv(1+Nk:end, :); % estimated neuropil
    F(:,irange)        = Fneu(:,irange) + Fdeconv(1:Nk, :); % estimated ROI signal
    
    ix = ix + NTbatch;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end

% get activity stats

sd = std(F, [], 2);

sk(:, 1) = skewness(F, [], 2);
sk(:, 2) = sd/mean(sd); 
sk(:, 3) = (max(F, [], 2)-median(F, 2))./sd;
sk(:, 4) = (prctile(F, 95, 2)-median(F, 2))./sd;

for j = 1:numel(stat)
   stat(j).dFstat = sk(j,:); 
end


end