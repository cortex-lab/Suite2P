function Fcell = get_signals_NEUmodel(ops, iplane)

try
   load(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat',...
       ops.ResultsSavePath,ops.mouse_name, ops.date, iplane, ops.Nk))
catch
   error('Could not find cell detection file \n') 
end

Nk = numel(stat);
Nkpar = ops.Nk;
%%
S = res.S;
M = res.M;
nBasis = round(sqrt(size(S,2)));

StS = S' * S;
LtL = zeros(Nk, Nk, 'single');
LtS = zeros(Nk, nBasis^2, 'single');
Ireg = diag([ones(Nk,1); zeros(nBasis^2,1)]);

for i = ops.Nk+1:Nk
%     ix = find(iclust==i);
    ix = stat(i).ipix;
    if numel(ix)==0
        LtL(i,i) = 0;
        LtS(i,:) = 0;
    else
        LtL(i,i) = sum(M(ix).^2);
        LtS(i,:) = M(ix) * S(ix, :);
    end
end

covL = [LtL LtS; LtS' StS];

covL = covL + 1e-4 * Ireg;
%% get signals  

[Ly Lx] = size(ops.mimg);

nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F        = zeros(Nk, sum(ops.Nframes), 'single');
Fneu    = zeros(Nk, sum(ops.Nframes), 'single');
% Fraw = zeros(Nk, sum(ops.Nframes), 'single');

covLinv = inv(covL((1+ops.Nk):end, (1+ops.Nk):end));
covLinv = covLinv ./ repmat(diag(covLinv), 1, size(covLinv,2));
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
       break; 
    end
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    NT= size(data,3);
    
    data = single(reshape(data, [], NT));
    %
%     data = data - repmat(mean(data,2), 1, size(data,2));
    
    StU = res.S' * data;
    %
    Ftemp = zeros(Nk, NT, 'single');
    for k = 1+ops.Nk:Nk
       ipix = stat(k).ipix; 
       if ~isempty(ipix)
           Ftemp(k,:) = res.M(ipix) * data(ipix,:);
       end
    end
    F(:,ix + (1:NT))        = Ftemp;
    
    i0 = (1+ops.Nk):Nk;
    %Fdeconv = covL((1+ops.Nk):end, (1+ops.Nk):end) \ cat(1, Ftemp(i0, :), StU);
    Fdeconv = covLinv * cat(1, Ftemp(i0, :), StU);
    
    Ftemp2 = Fdeconv(1:(Nk-ops.Nk),:);
    
    Fneu(i0,ix + (1:NT))     = Ftemp2;
    
        
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);
%%
csumNframes = [0 cumsum(ops.Nframes)];
Fcell       = cell(1, length(ops.Nframes));
FcellNeu    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i} = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeu{i} = Fcell{i} - F(:, csumNframes(i) + (1:ops.Nframes(i)));
%     FcellNeu{i} =  F(:, csumNframes(i) + (1:ops.Nframes(i))) - Fcell{i};
end

save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, ops.Nk),  'ops', 'res', 'stat', ...
    'stat0', 'res0', 'Fcell', 'FcellNeu')

