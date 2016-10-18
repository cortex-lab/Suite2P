function mimgR = red_channel_mean(ops)

% numPlanes = length(ops.planesToProcess);

%% build file list with red channel

if (isfield(ops, 'SubDirsRed') && ~isempty(ops.SubDirsRed))
   subDirsRed = ops.SubDirsRed;
else
    if (isfield(ops, 'expred') && ~isempty(ops.expred))
        for i = 1:length(ops.expred)
            subDirsRed{i} = sprintf('%d', ops.expred(i));
        end
    else
        warning('could not find red channel info, returning...')
        return;
    end 
end
%%
% build file list
fsrootRED = [];

for j = 1:length(subDirsRed)
    fsrootRED{j} = dir(fullfile(ops.RootDir, subDirsRed{j}, '*.tif'));
    jn=find(ops.expts==str2num(subDirsRed{j}));
    if ~isempty(jn)
        for k = 1:length(fsrootRED{jn})
            fsrootRED{j}(k).name = fullfile(ops.RootDir, subDirsRed{j}, fsrootRED{j}(k).name);
        end
    end
end

%% check if the RED channel file has already been registered, find indices to shift by
fsRED   = [];
fs      = [];
nimg    = [];
nimgall = [];

nimgFirst = zeros(length(ops.fsroot), ops.nplanes);

for j = 1:length(ops.fsroot)
    if ~isempty(ops.fsroot{j})
        ops.fsroot{j}(1).nFrames = nFrames(ops.fsroot{j}(1).name);
        
        for k = 1:length(ops.fsroot{j})
            if abs(ops.fsroot{j}(1).bytes - ops.fsroot{j}(k).bytes)<10
                ops.fsroot{j}(k).nFrames = ops.fsroot{j}(1).nFrames;
            else
                ops.fsroot{j}(k).nFrames = nFrames(ops.fsroot{j}(k).name);
            end
        end
        
        fs = [fs {ops.fsroot{j}.name}];
        nimg{j} = [ops.fsroot{j}.nFrames];
        nimgall = [nimgall ops.fsroot{j}.nFrames];
    end
end
%
nindx = zeros(1, ops.nplanes);
for k = 1:length(fs)
    nimgFirst(k, :) = nindx;
    nindx = nindx + ceil((-[0:1:(ops.nplanes-1)] + nimgall(k)/ops.nchannels_red)/ops.nplanes);
end
%%
for j = 1:length(fsrootRED)
    fsRED = [fsRED {fsrootRED{j}.name}];
end
indx = zeros(length(fsRED), 1);
for k = 1:length(fsRED)
    indx(k) = find(cellfun(@(x) strcmp(x, fsRED{k}), fs));
end

% sum(nimgall(indx))/ops.nchannels/ops.nplanes
ntifs = ceil(500/length(fsRED));

%
DS = cell(ops.nplanes, 1);
for j = 1:ops.nplanes
   try 
       root = ops.ResultsSavePath;
       fname = sprintf('regops_%s_%s_plane%d.mat', ops.mouse_name, ops.date, j);
       load(fullfile(root, fname))
       DS{j} = ops.DS;
   catch       
   end
end
%
ntf0 = 0;
numPlanes = ops.nplanes;
iplane0 = 1:1:ops.nplanes;
for k = 1:length(fsRED)
    if nimgall(indx(k))>=median(nimgall(indx))        
        iplane0 = mod(iplane0-1, numPlanes) + 1;
         
        ichanset = [ops.nchannels_red*ops.nplanes + [ops.nchannels_red (ntifs*ops.nchannels_red*ops.nplanes)] ...
            ops.nchannels_red]; 
        data = loadFrames(fsRED{k}, ichanset(1),ichanset(2), ichanset(3));        
        if ~exist('mimgR')
            [Ly, Lx, ~] = size(data);
            mimgR = zeros(Ly, Lx, ops.nplanes);
        end
        data = reshape(data, Ly, Lx, ops.nplanes, ntifs);         

        for j = 1:size(data,3)
            dsall = DS{j}(nimgFirst(indx(k), j)+1 + [1:size(data,4)], :);
            data(:, :, j,:)        = ...
                    register_movie(data(:, :, j, :), ops, dsall);
        end
        mimgR = mimgR + mean(data(:,:,iplane0,:), 4);            
        ntf0 = ntf0 + 1;
        
        nFr = img.nFrames(fsRED{k});
        iplane0 = iplane0 - nFr/ops.nchannels_red;
        
%         keyboard;
    end
end

mimgR = mimgR/ntf0;