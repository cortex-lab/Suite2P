function ops1 = red_channel_mean(ops)

numPlanes = length(ops.planesToProcess);

%% build file list with red channel

% build file list
ops.fsrootRED = [];

for j = 1:length(ops.SubDirsRed)
    ops.fsrootRED{j} = dir(fullfile(ops.RootDir, ops.SubDirsRed{j}, '*.tif'));
    for k = 1:length(ops.fsroot{j})
        ops.fsrootRED{j}(k).name = fullfile(ops.RootDir, ops.SubDirsRed{j}, ops.fsrootRED{j}(k).name);
    end
end

fs = ops.fsrootRED;
ntifs = sum(cellfun(@(x) numel(x), fs));

nfmax = ceil(2000/ntifs);

%% check if the RED channel file has already been registered, find indices to shift by


%% find plane