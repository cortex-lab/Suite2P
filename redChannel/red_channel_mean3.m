function [mimgR,mimgG] = red_channel_mean3(ops)

% numPlanes = length(ops.planesToProcess);

% build file list with red channel

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
%
% build file list
fsRED = dir(fullfile(ops.RootDir, subDirsRed{1}, '*.tif'));
for k = 1:length(fsRED)
    fsRED(k).name = fullfile(ops.RootDir, subDirsRed{1}, fsRED(k).name);
end

ops1 = cell(ops.nplanes, 1);
for j = 1:ops.nplanes
    root = ops.ResultsSavePath;
    fname = sprintf('regops_%s_%s_plane%d.mat', ops.mouse_name, ops.date, j);
    dat = load(fullfile(root, fname));
    ops1{j} = dat.ops;
    ops1{j}.useGPU = ops.useGPU;
end

%
ntf0 = 0;
numPlanes = ops.nplanes;
iplane0 = 1:1:ops.nplanes;
for k = 1:length(fsRED)
    iplane0 = mod(iplane0-1, numPlanes) + 1;
    nFr = nFramesTiff(fsRED(k).name);
    data = loadFramesBuff(fsRED(k).name, 1, nFr, 1, ops.temp_tiff);
    
    if ~exist('mimgR', 'var')
        [Ly, Lx, ~] = size(data);
        mimgR = zeros(Ly, Lx, ops.nplanes);
    end
    if ~exist('mimgG', 'var')
        [Ly, Lx, ~] = size(data);
        mimgG = zeros(Ly, Lx, ops.nplanes);
    end
    %
    data = data(:,:,1:floor(nFr/(2*ops.nplanes))*2*ops.nplanes);
    data = reshape(data, Ly, Lx, 2,ops.nplanes, []);
    dataG = sq(data(:,:,1,:,:), 3);
    dataR = sq(data(:,:,2,:,:), 3);
    
    for i = 1:numPlanes
        [ds, ~]  = registration_offsets(squeeze(dataG(:,:,i,:)), ops1{i}, 0);
        if k==1
            ds(1,:) = 0;
        end
        dataR(:, :, i,:)        = ...
            register_movie(squeeze(dataR(:, :, i, :)), ops1{i}, ds);
        dataG(:, :, i,:)        = ...
            register_movie(squeeze(dataG(:, :, i, :)), ops1{i}, ds);
    end
    
    mimgR = mimgR + mean(dataR(:,:,iplane0,:), 4);
    mimgG = mimgG + mean(dataG(:,:,iplane0,:), 4);
    ntf0 = ntf0 + 1;
    
    iplane0 = iplane0 - nFr/ops.nchannels_red;
end

mimgR = mimgR/ntf0;
