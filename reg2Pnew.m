function ops1 = reg2Pnew(ops)
%%
% iplane = ops.iplane;
% bitspersamp = 16;
numPlanes = length(ops.planesToProcess);

if isfield(ops, 'chunk_align') && ~isempty(ops.chunk_align); chunk_align   = ops.chunk_align(1);
else chunk_align = 1; end
if isfield(ops, 'nplanes') && ~isempty(ops.nplanes); nplanes   = ops.nplanes;
else nplanes = 1; end
if isfield(ops, 'nchannels') && ~isempty(ops.nchannels); nchannels   = ops.nchannels;
else nchannels = 1; end
if isfield(ops, 'gchannel') && ~isempty(ops.gchannel); ichannel    = ops.gchannel;
else ichannel = 1; end
if nchannels>2 && isfield(ops, 'rchannel') && ~isempty(ops.rchannel); rchannel   = ops.rchannel;
else rchannel = 2; end
if isfield(ops, 'AlignToRedChannel') && ~isempty(ops.AlignToRedChannel); red_align = ops.AlignToRedChannel;
else red_align = 0; end
if isfield(ops, 'BiDiPhase') && ~isempty(ops.BiDiPhase); BiDiPhase = ops.BiDiPhase;
else BiDiPhase = 0; end
if isfield(ops, 'LoadRegMean') && ~isempty(ops.LoadRegMean); LoadRegMean = ops.LoadRegMean;
else LoadRegMean = 0; end

fs = ops.fs;

%% find the mean frame after aligning a random subset
ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = floor(ops.NimgFirstRegistration/ntifs);
if nfmax>=2000
    nfmax = 1999;
end
nbytes = fs{1}(1).bytes;
nFr = nFrames(fs{1}(1).name);
Info0 = imfinfo(fs{1}(1).name);
Ly = Info0(1).Height;
Lx = Info0(1).Width;
ops.Lx = Lx;
ops.Ly = Ly;

indx = 0;
IMG = zeros(Ly, Lx, nplanes, ops.NimgFirstRegistration, 'single');

if ops.doRegistration
    for k = 1:length(ops.SubDirs)
        iplane0 = 1;
        for j = 1:length(fs{k})
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = nFrames(fs{k}(j).name);
            end
            if nFr<(nchannels*nplanes*nfmax + nchannels*nplanes)
                continue;
            end
            
            iplane0 = mod(iplane0-1, nplanes) + 1;
            if red_align
                ichanset = [nchannels*nplanes + nchannels*(iplane0-1) + [rchannel;...
                    nchannels*nplanes*nfmax]; nchannels];
            else
                ichanset = [nchannels*nplanes + nchannels*(iplane0-1) + [ichannel;...
                    nchannels*nplanes*nfmax]; nchannels];
            end
            iplane0 = iplane0 - nFr/nchannels;
            data = loadFrames(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
            data = reshape(data, Ly, Lx, nplanes, []);
            
            if BiDiPhase
                yrange = 2:2:Ly;
                if BiDiPhase>0
                    data(yrange, (1+BiDiPhase):Lx,:,:) = data(yrange, 1:(Lx-BiDiPhase),:,:);
                else
                    data(yrange, 1:Lx+BiDiPhase,:,:) = data(yrange, 1-BiDiPhase:Lx,:,:);
                end
            end
            IMG(:,:,:,indx+(1:size(data,4))) = data;
            indx = indx + size(data,4);
            
        end
    end
    IMG =  IMG(:,:,:,1:indx);
    %
    for i = 1:numPlanes
        ops1{i} = align_iterative(squeeze(IMG(:,:,ops.planesToProcess(i),:)), ops);
    end
end

if ops.showTargetRegistration
    figure('position', [900 50 900 900])
    ax = ceil(sqrt(numel(ops1)));
    for i = 1:length(ops1)
        subplot(ax,ax,i)
        imagesc(ops1{i}.mimg)
        colormap('gray')
        title(sprintf('Registration for plane %d, mouse %s, date %s', ...
            i, ops.mouse_name, ops.date))
    end
    
    drawnow
end
%
clear IMG

%%
for i = 1:numPlanes
    ops1{i}.RegFile = fullfile(ops.RegFileRoot, sprintf('tempreg_plane%d.bin', ops.planesToProcess(i)));
    regdir = fileparts(ops1{i}.RegFile);
    if ~exist(regdir, 'dir')
        mkdir(regdir);
    end    
    
    % open bin file for writing
    fid{i}             = fopen(ops1{i}.RegFile, 'w');
    ops1{i}.DS          = [];
    ops1{i}.CorrFrame   = [];
    ops1{i}.mimg1       = zeros(Ly, Lx);
    
end
nbytes = fs{1}(1).bytes;
nFr = nFrames(fs{1}(1).name);
if ops.useImRead
    Info = imfinfo(fs{1}(1).name);
end

tic
for k = 1:length(fs)
    for i = 1:numPlanes
         ops1{i}.Nframes(k)  = 0;
    end
    
    iplane0 = 1:1:ops.nplanes;
    for j = 1:length(fs{k})
        iplane0 = mod(iplane0-1, numPlanes) + 1;
        if ops.useImRead
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = nFrames(fs{k}(j).name);
                Info = imfinfo(fs{k}(j).name);
            end
            if red_align
%                 ichanrange = (nchannels*(iplane0-1)+rchannel):nchannels:nFr;
                ichanrange = rchannel:nchannels:nFr;
            else
                ichanrange = ichannel:nchannels:nFr;
%                 ichanrange = (nchannels*(iplane0-1)+ichannel):nchannels:nFr;
            end
            data = zeros(Ly, Lx, numel(ichanrange), ops.RawPrecision);
            for ix = 1:length(ichanrange)
                data(:,:,ix) = imread(fs{k}(j).name, 'Index', ichanrange(ix),'Info', Info);
            end
        else
            nFr = nFrames(fs{k}(j).name);
            if red_align
%                 ichanset = [nchannels*(iplane0-1)+rchannel; nFr; nchannels];
                ichanset = [rchannel; nFr; nchannels];
            else
                ichanset = [ichannel; nFr; nchannels];
            end
            data = loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
        end
        
        if BiDiPhase
            yrange = 2:2:Ly;
            if BiDiPhase>0
                data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
            else
                data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
            end
        end
        
        if ops.doRegistration
            % get the registration offsets
            dsall = zeros(size(data,3), 2);
            for i = 1:numPlanes
                ifr0 = iplane0(ops.planesToProcess(i));
                indframes = ifr0:nplanes:size(data,3);
                [ds, Corr]  = registration_offsets(data(:,:,indframes), ops1{i}, 0);
                dsall(indframes,:)  = ds;
                % collect ds
                if j==1
                    ds(1,:) = 0;
                end
                ops1{i}.DS          = cat(1, ops1{i}.DS, ds);
                ops1{i}.CorrFrame   = cat(1, ops1{i}.CorrFrame, Corr);
            end
            
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if red_align
                if ops.useImRead
%                     ichanrange = (nchannels*(iplane0-1)+ichannel):nchannels:length(Info);
                    ichanrange = ichannel:nchannels:length(Info);
                    data = zeros(Ly, Lx, numel(ichanrange), ops.RawPrecision);
                    for ix = 1:length(ichanrange)
                        data(:,:,ix) = imread(fs{k}(j).name, 'Index', ichanrange(ix),'Info', Info);
                    end
                else
                    ichanset = [ichannel; nFr; nchannels];
                    data = loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
                end
            end
            
            ix0 = 0;
            Nbatch = 1000;
            dreg = zeros(size(data), ops.RegPrecision);
            while ix0<size(data,3)
                indxr = ix0 + (1:Nbatch);
                indxr(indxr>size(data,3)) = [];
                dreg(:, :, indxr)        = ...
                    register_movie(data(:, :, indxr), ops1{1}, dsall(indxr,:));
                ix0 = ix0 + Nbatch;
            end
           
        else
            dreg = data;
        end
        % write dreg to bin file
        for i = 1:numPlanes
            ifr0 = iplane0(ops.planesToProcess(i));
            indframes = ifr0:nplanes:size(data,3);
            dwrite = dreg(:,:,indframes);
            fwrite(fid{i}, dwrite, ops.RegPrecision);
            
            ops1{i}.Nframes(k) = ops1{i}.Nframes(k) + size(dwrite,3);
            ops1{i}.mimg1 = ops1{i}.mimg1 + sum(dwrite,3);
        end
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        % delete temporarily copied tiffs
        if ops.CopyDataLocally && ops.DeleteRawOnline
            % check if the location is NOT on zserver
            if ~isempty(strfind(ops.TempStorage, 'zserver')) || ...
                    strcmp(ops.TempStorage(1), '\') || ...
                    strcmp(ops.TempStorage(1), '/')
                warning('You are trying to remove a file from a network location, skipping...')
            else
                delete(fs{k}(j).name);
            end
        end
        
        iplane0 = iplane0 - nFr/nchannels;
    end
    for i = 1:numPlanes
        ops1{i}.mimg1 = ops1{i}.mimg1/ops1{i}.Nframes(k);
    end
    
end
% close bin file
%%
for i = 1:numPlanes    
    fclose(fid{i});
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i);
    end
    
end
%%
% compute xrange, yrange
for i = 1:numPlanes
    if ops.doRegistration
        minDs = min(ops1{i}.DS(2:end, [1 2]), [], 1);
        maxDs = max(ops1{i}.DS(2:end, [1 2]), [], 1);
        disp([minDs(1) maxDs(1) minDs(2) maxDs(2)])
        if BiDiPhase>0
            maxDs(2) = max(1+BiDiPhase, maxDs(2));
        elseif BiDiPhase<0
            minDs(2) = min(BiDiPhase, minDs(2));
        end
        
        ops1{i}.yrange = ceil(maxDs(1)):floor(Ly+minDs(1));
        ops1{i}.xrange = ceil(maxDs(2)):floor(Lx+minDs(2));
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end
    
    savepath = sprintf('%s/', ops.ResultsSavePath);
    
    if ~exist(savepath, 'dir')
        mkdir(savepath)
    end
    ops = ops1{i};
    save(sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date,  ops.planesToProcess(i)),  'ops')
end

% delete temporarily copied tiffs
if ops.CopyDataLocally
    % check if the location is NOT on zserver
    if ~isempty(strfind(ops.TempStorage, 'zserver')) || ...
            strcmp(ops.TempStorage(1), '\') || ...
            strcmp(ops.TempStorage(1), '/')
        warning('You are trying to remove a file from a network location, skipping...')
    else
        rmdir(fullfile(ops.TempStorage, ops.mouse_name), 's');
        %             rmdir(fullfile(ops.TempDir), 's');
    end
end

%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
