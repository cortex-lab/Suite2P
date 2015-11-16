function [R, ops] = PixelMapsSVD(ton, block, SVDfilepath, ops1)
% ton is a cell, each containing stimulus start times for many repetitions
% of the same stimulus

load(SVDfilepath)


% this adds the options from ops1 to ops
ops.ExpFiltTau      = ops1.ExpFiltTau;
ops.ntf             = ops1.ntf; % timepoints for the filter
ops.method          = ops1.method; % deconv or sta
ops.nSVDforPixelMaps = ops1.nSVDforPixelMaps;
%%
% root = sprintf('%s//%s//%s', ops.ResultsSavePath, mouse_name, sessid);
% fname = sprintf('SVD_%s_%s_plane%d.mat', ops.mouse_name,ops.date, ops.iplane);
% load(fullfile(root, fname));

%%

iblock = strcmp(ops.SubDirs, num2str(block));
F = Vcell{iblock};

% F at this point should be nSVD by NT
F = F - expfilt(F, ops.ExpFiltTau);

[R, filt, STA] = get_responses3(F, ton, ops);

%%
nSVD = ops.nSVDforPixelMaps;
pixR = reshape(U(:, :, 1:nSVD), [], nSVD) * STA(1:nSVD, :);

pixR = reshape(pixR, size(U,1), size(U,2), size(STA,2), size(STA,3));

RR = reshape(pixR, [], size(pixR,4));
RR = RR'*RR;

[filt, Sv]=eigs(double(RR), 1);
[~, imaxfilt] = max(abs(filt));
filt = filt * sign(filt(imaxfilt));

R = reshape(pixR, [], size(pixR,4)) * filt;
R = reshape(R, size(U,1), size(U,2), []);


%%