function [ops, stat, Fcell, FcellNeu]      = extractSignals(ops, m, stat)

Nk = numel(stat);

S = reshape(m.S, [], size(m.S, ndims(m.S)));

nBasis = size(S,2);

Ireg = diag([ones(Nk,1); zeros(nBasis,1)]);

% maskNeu = ones(size(S,1), 1);

covL = [m.LtL m.LtS; m.LtS' m.StS];

covLinv = inv(covL + 1e-4 * Ireg);
covLinv = covLinv ./ repmat(diag(covLinv), 1, size(covLinv,2));
%% get signals  

nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F        = zeros(Nk, sum(ops.Nframes), 'single');
Fneu    = zeros(Nk, sum(ops.Nframes), 'single');

% S    = bsxfun(@times, S, maskNeu(:));

[Ly Lx] = size(ops.mimg1);

% find indices of good clusters 
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
       break; 
    end
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    NT= size(data,3);
    
    % process the data
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
    Fdeconv = covLinv * cat(1, Ftemp, S' * data);
    
    Fneu(:,ix + (1:NT))     = m.LtS * Fdeconv(1+Nk:end, :); % estimated neuropil
    F(:,ix + (1:NT))        = Fneu(:,ix + (1:NT)) + Fdeconv(1:Nk, :); % estimated ROI signal
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
    
end
fclose(fid);
%% get activity stats

sd = std(F, [], 2);

sk(:, 1) = skewness(F, [], 2);
sk(:, 2) = sd/mean(sd); 
sk(:, 3) = (max(F, [], 2)-median(F, 2))./sd;
sk(:, 4) = (prctile(F, 95, 2)-median(F, 2))./sd;

for j = 1:numel(stat)
   stat(j).dFstat = sk(j,:); 
end

%%
csumNframes = [0 cumsum(ops.Nframes)];
Fcell       = cell(1, length(ops.Nframes));
FcellNeu    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i}     = F(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeu{i}  = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end

% save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
%     ops.mouse_name, ops.date, ops.iplane),  'ops',  'stat',...
%         'Fcell', 'FcellNeu', 'clustrules', '-v7.3')


