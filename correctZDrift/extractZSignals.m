function [FcellZ, FcellNeuZ] = extractZSignals(ops, zpos, Fcell, FcellNeu, subpixel)

Nk = size(Fcell{1},1);

F = [];
Fneu = [];
for k = 1:numel(Fcell)
    F = [F Fcell{k}];
    Fneu = [Fneu FcellNeu{k}];
end

% compute each cell's profile in Z
zspread  = 10;
iZ = [-zspread : 1/subpixel : zspread];
zind = zpos * subpixel + find(iZ==0);
zF = zeros(Nk, numel(iZ));
zFneu = zeros(Nk, numel(iZ));
for j = 1:numel(iZ)
    zF(:,j)    = mean(F(:, zind == j), 2);
    zFneu(:,j) = mean(Fneu(:, zind == j), 2);
end

if sum(zind==1) == 0
    jm = min(zind);
    zF(:, 1:jm-1) = repmat(zF(:, jm), 1, jm-1);
    zFneu(:, 1:jm-1) = repmat(zFneu(:, jm), 1, jm-1);
end

if sum(zind==numel(iZ)) == 0
    jm = max(zind);
    zF(:, jm+1:end) = repmat(zF(:, jm), 1, numel(iZ) - jm);
    zFneu(:, jm+1:end) = repmat(zFneu(:, jm), 1, numel(iZ) - jm);
end

zF= fixnangaps(zF);
zFneu= fixnangaps(zFneu);

zF = my_conv2(zF, 0.5, 2);
zFneu = my_conv2(zFneu, 0.5, 2);
plot(zF(1,:));

%%
Fz      = zF(:, zind);
Fzneu    = zFneu(:, zind);    
    
csumNframes = [0 cumsum(ops.Nframes)];
FcellZ       = cell(1, length(ops.Nframes));
FcellNeuZ    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    FcellZ{i}     = Fz(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeuZ{i}  = Fzneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end