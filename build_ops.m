function ops = build_ops(db, ops)

% ops = db;
fieldNames = fieldnames(db);
for i = 1:size(fieldNames,1)
    ops.(fieldNames{i}) = db.(fieldNames{i});
end

if isfield(ops, 'chunk_align') && ~isempty(ops.chunk_align)
    if numel(ops.chunk_align)<ops.nplanes
        if numel(ops.chunk_align)==1
            ops.chunk_align = repmat(ops.chunk_align, 1, ops.nplanes);
        else
            warning('Incompatible number of inputs in chunk_align, filling up with ones \n')
            ops.chunk_align(end+1:ops.nplanes) = 1;
        end
    end
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
end

% ops.nchannels  = db.nchannels;
% ops.nplanes    = db.nplanes;
% if ops.nchannels>1
%     ops.gchannel   = db.gchannel;
% else
%     ops.gchannel   = 1;
% end

% ops.mouse_name = db.mouse_name;
% ops.date       = db.date ;


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
if ops.CopyDataLocally
    ops.fs = [];
    for j = 1:length(ops.SubDirs)
        ops.fs{j} = dir(fullfile(ops.TempDir, ops.SubDirs{j}, '*.tif'));
        for k = 1:length(ops.fs{j})
            ops.fs{j}(k).name = fullfile(ops.TempDir, ops.SubDirs{j}, ops.fs{j}(k).name);
        end
    end
else
    ops.fs = ops.fsroot;
end