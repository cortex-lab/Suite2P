function ops = build_ops3(db, ops)

% ops = db;
ops = addfields(ops, db);

for k = 1:length(db.expts)
    ops.SubDirs{k}    = num2str(db.expts(k));
end

ops.RootDir = fullfile(ops.RootStorage, ops.mouse_name, ops.date);

% build file list
ops.fsroot = [];
for j = 1:length(ops.SubDirs)
    ops.fsroot{j} = dir(fullfile(ops.RootDir, ops.SubDirs{j}, '*.tif'));
    for k = 1:length(ops.fsroot{j})
        ops.fsroot{j}(k).name = fullfile(ops.RootDir, ops.SubDirs{j}, ops.fsroot{j}(k).name);
    end
end

try    
     % MK code for automatically determining number of planes and channels
    [~, header] = loadFramesBuff(ops.fsroot{1}(1).name, 1, 1, 1);
    
    hh=header{1};
    str = hh(strfind(hh, 'channelsSave = '):end);
    ind = strfind(str, 'scanimage');
    ch = str2num(str(16 : ind(1)-1));
    ops.nchannels = length(ch);
    
    fastZEnable = sscanf(hh(findstr(hh, 'fastZEnable = '):end), 'fastZEnable = %d');
    fastZDiscardFlybackFrames = sscanf(hh(findstr(hh, 'fastZDiscardFlybackFrames = '):end), 'fastZDiscardFlybackFrames = %d');
    if isempty(fastZDiscardFlybackFrames)
        fastZDiscardFlybackFrames = 0;
    end
    stackNumSlices = sscanf(hh(findstr(hh, 'stackNumSlices = '):end), 'stackNumSlices = %d');
    
    ops.nplanes = 1;
    if fastZEnable
        ops.nplanes = stackNumSlices+fastZDiscardFlybackFrames;
    end
    
    str = hh(strfind(hh, 'scanZoomFactor = '):end);
    ind = strfind(str, 'scanimage');
    ops.zoomMicro = str2double(str(18 : ind(1)-1));
catch
end

if ~(isfield(ops, 'planesToProcess') && ~isempty(ops.planesToProcess))
    ops.planesToProcess = 1:ops.nplanes;
else
    % planesToProcess is not working right now
    ops.planesToProcess = 1:ops.nplanes; 
end

CharSubDirs = '';
for i = 1:length(ops.SubDirs)
    CharSubDirs = [CharSubDirs ops.SubDirs{i} '_'];
end
CharSubDirs = CharSubDirs(1:end-1);

ops.ResultsSavePath = sprintf('%s//%s//%s//%s//', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
        CharSubDirs);