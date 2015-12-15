function ops = build_ops2(db, ops)

% ops = db;
fieldNames = fieldnames(db);
for i = 1:size(fieldNames,1)
    if ~isempty(db.(fieldNames{i}))
        ops.(fieldNames{i}) = db.(fieldNames{i});
    end
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
else
    % planesToProcess is not working right now
    ops.planesToProcess = 1:ops.nplanes; 
end

for k = 1:length(db.expts)
    ops.SubDirs{k}    = num2str(db.expts(k));
end

ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);
ops.TempDir = fullfile(ops.TempStorage, ops.mouse_name, ops.date);


% build file list
ops.fsroot = [];
for j = 1:length(ops.SubDirs)
    ops.fsroot{j} = dir(fullfile(ops.RootDir, ops.SubDirs{j}, '*.tif'));
    for k = 1:length(ops.fsroot{j})
        ops.fsroot{j}(k).name = fullfile(ops.RootDir, ops.SubDirs{j}, ops.fsroot{j}(k).name);
    end
end

% build file list
ops.fs = [];
for j = 1:length(ops.SubDirs)
    ops.fs{j} = dir(fullfile(ops.TempDir, ops.SubDirs{j}, '*.tif'));
    for k = 1:length(ops.fs{j})
        if ops.CopyDataLocally
            ops.fs{j}(k).name = fullfile(ops.TempDir, ops.SubDirs{j}, ops.fs{j}(k).name);
        else
            ops.fs{j}(k).name =  ops.fsroot{j}(k).name;
        end
    end
end

CharSubDirs = '';
for i = 1:length(ops.SubDirs)
    CharSubDirs = [CharSubDirs ops.SubDirs{i} '_'];
end
CharSubDirs = CharSubDirs(1:end-1);

ops.ResultsSavePath = sprintf('%s//%s//%s//%s//', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
        CharSubDirs);