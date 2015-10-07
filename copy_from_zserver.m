function copy_from_zserver(ops)

if ~exist(ops.TempDir, 'dir');  mkdir(ops.TempDir); end
for k = 1:length(ops.SubDirs)
    if ~exist(fullfile(ops.TempDir, ops.SubDirs{k}), 'dir')
        copyfile(fullfile(ops.RootDir, ops.SubDirs{k}), fullfile(ops.TempDir, ops.SubDirs{k}))
    end
end
