function ops = reg2P(ops)
%%
iplane = ops.iplane;
bitspersamp = 16;

if isfield(ops, 'chunk_align') && ~isempty(ops.chunk_align); chunk_align   = ops.chunk_align(iplane);
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

addpath(genpath('\\zserver\Code\Register\'))
fs = ops.fs;

%% find the mean frame after aligning a random subset

ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = chunk_align*floor(ops.NimgFirstRegistration/ntifs);
if nfmax>=2000
    nfmax = 2000-chunk_align;
end

nbytes = fs{1}(1).bytes;
nFr = img.nFrames(fs{1}(1).name);
Info0 = imfinfo(fs{1}(1).name);

Ly = Info0(1).Height;
Lx = Info0(1).Width;

if ops.doRegistration
    flag = 1;
    if LoadRegMean || (isfield(ops, 'regopspath') && ~isempty(ops.regopspath))
        flag = 0;
        if ~isfield(ops, 'regopspath')
            fname = sprintf('%s/%s/%s/regops_%s_%s_plane%d.mat', ...
                ops.ResultsSavePath, ops.mouse_name, ops.date, ...
                ops.mouse_name, ops.date, iplane);
        else
            fname = sprintf('%s%d.mat', ops.regopspath, iplane);
        end
        if ~exist(fname, 'file')
            warning('Options file does not exist at %s\n recomputing target frame...', fname);
            flag = 1;
        else
            dd = load(fname);
            ops.mimg 	= dd.ops.mimg;
        end
    end
    if flag
        IMG = zeros(Ly, Lx, ops.NimgFirstRegistration, 'single');
        
        indx = 0;
        for k = 1:length(ops.SubDirs)
            iplane0 = iplane;
            for j = 1:length(fs{k})
                if abs(nbytes - fs{k}(j).bytes)>1e3
                    nbytes = fs{k}(j).bytes;
                    nFr = img.nFrames(fs{k}(j).name);
                end
                if nFr<(nchannels*nplanes*nfmax + nchannels*nplanes)
                    continue;
                end
                
                iplane0 = mod(iplane0-1, nplanes) + 1;
                if red_align
                    ichanset = [nchannels*nplanes + [nchannels*(iplane0-1)+rchannel;...
                        nchannels*nplanes*nfmax]; nplanes*nchannels];
                else
                    ichanset = [nchannels*nplanes + [nchannels*(iplane0-1)+ichannel;...
                        nchannels*nplanes*nfmax]; nplanes*nchannels];
                end
                iplane0 = iplane0 - nFr/nchannels;
                data = img.loadFrames(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
                data = squeeze(mean(reshape(data, Ly, Lx, chunk_align, []), 3));
                
                if BiDiPhase
                    yrange = 2:2:Ly;
                    if BiDiPhase>0
                        data(yrange, (1+BiDiPhase):Lx,:) = data(yrange, 1:(Lx-BiDiPhase),:);
                    else
                        data(yrange, 1:Lx+BiDiPhase,:) = data(yrange, 1-BiDiPhase:Lx,:);
                    end
                end
                IMG(:,:,indx+(1:size(data,3))) = data;
                indx = indx + size(data,3);
                
            end
        end
        IMG =  IMG(:,:,1:indx);
        %
        ops = align_iterative(IMG, ops);
    end
		
    if ops.showTargetRegistration
        figure;
        imagesc(ops.mimg)
        colormap('gray')
        title(sprintf('Registration for plane %d, mouse %s, date %s', iplane, ops.mouse_name, ops.date))
        drawnow
    end
    %
    clear IMG
end
ops.RegFile = [ops.RegFileRoot, sprintf('_plane%d.bin', iplane)];
regdir = fileparts(ops.RegFile);
if ~exist(regdir, 'dir')
   mkdir(regdir); 
end


% keyboard;

%%
ops.DS = [];
% open bin file for writing
fid             = fopen(ops.RegFile, 'w');
ops.DS          = [];
ops.CorrFrame   = [];
ops.mimg1       = zeros(Ly, Lx);
nbytes = fs{1}(1).bytes;
nFr = img.nFrames(fs{1}(1).name);
if ops.useImRead
   Info = imfinfo(fs{1}(1).name); 
end

tic
for k = 1:length(fs)
    ops.Nframes(k)  = 0;
    iplane0 = iplane;
    for j = 1:length(fs{k})
        iplane0 = mod(iplane0-1, nplanes) + 1;
        if ops.useImRead
            if abs(nbytes - fs{k}(j).bytes)>1e3
                nbytes = fs{k}(j).bytes;
                nFr = img.nFrames(fs{k}(j).name);
                Info = imfinfo(fs{k}(j).name);
            end
            if red_align
                ichanrange = (nchannels*(iplane0-1)+rchannel):nplanes*nchannels:nFr;
            else
                ichanrange = (nchannels*(iplane0-1)+ichannel):nplanes*nchannels:nFr;
            end
            data = zeros(Ly, Lx, numel(ichanrange), ops.RawPrecision);
            for ix = 1:length(ichanrange)
                data(:,:,ix) = imread(fs{k}(j).name, 'Index', ichanrange(ix),'Info', Info);
            end
        else
            nFr = img.nFrames(fs{k}(j).name);
            if red_align
                ichanset = [nchannels*(iplane0-1)+rchannel; nFr; nplanes*nchannels];
            else
                ichanset = [nchannels*(iplane0-1)+ichannel; nFr; nplanes*nchannels];
            end
            data = img.loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
        end
        iplane0 = iplane0 - nFr/nchannels;
        
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
            NTb = size(data,3);
            NTb0 =  chunk_align * floor(NTb/chunk_align);
            data0 = single(data(:, :,1:NTb0));
            data0 = squeeze(mean(reshape(data0, Ly, Lx, chunk_align, []), 3));
            [ds, Corr]  = registration_offsets(data0, ops, 0);
            ds = reshape(repmat(ds', chunk_align, 1), 2, [])';
            ds(NTb0+1:NTb, :) = repmat(ds(NTb0, :), NTb - NTb0, 1);
            
            % if aligning by the red channel, data needs to be reloaded as the
            % green channel
            if red_align
                if ops.useImRead
                    ichanrange = (nchannels*(iplane-1)+ichannel):nplanes*nchannels:length(Info);
                    data = zeros(Ly, Lx, numel(ichanrange), ops.RawPrecision);
                    for ix = 1:length(ichanrange)
                        data(:,:,ix) = imread(fs{k}(j).name, 'Index', ichanrange(ix),'Info', Info);
                    end
                else
                    ichanset = [nchannels*(iplane-1)+ichannel; nFr; nplanes*nchannels];
                    data = img.loadFrames(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
                end
            end
            
            ix0 = 0;
            Nbatch = 1000;
            dreg = zeros(size(data), ops.RegPrecision);
            while ix0<size(data,3)
                indxr = ix0 + (1:Nbatch);
                indxr(indxr>size(data,3)) = [];
                dreg(:, :, indxr)        = ...
                    register_movie(data(:, :, indxr), ops, ds(indxr,:));
                ix0 = ix0 + Nbatch;
            end
            % collect ds
            ops.DS          = cat(1, ops.DS, ds);
            ops.CorrFrame   = cat(1, ops.CorrFrame, Corr);
        else
            dreg = data;
        end
        % write dreg to bin file
        fwrite(fid, dreg, ops.RegPrecision);        
        
        ops.Nframes(k) = ops.Nframes(k) + size(data,3);
        ops.mimg1 = ops.mimg1 + sum(dreg,3);        
        
        if ~isempty(ops.RegFileTiffLocation)
            [~, partname , ~] = fileparts(ops.fs{k}(j).name);           
            foldr = fullfile(ops.RegFileTiffLocation, ops.mouse_name, sprintf('Plane%d', iplane), ...
                ops.SubDirs{k});
            if ~exist(foldr, 'dir')
                mkdir(foldr)
            end
            fname = fullfile(foldr, partname);            
            TiffWriter(uint16(dreg),[fname '.tif'],bitspersamp); 
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
    end
    ops.mimg1 = ops.mimg1/ops.Nframes(k);
    
end
% close bin file
fclose(fid);

% compute xrange, yrange
if ops.doRegistration
    minDs = min(ops.DS(2:end, [1 2]), [], 1);
    maxDs = max(ops.DS(2:end, [1 2]), [], 1);
    if BiDiPhase>0
        maxDs(2) = max(1+BiDiPhase, maxDs(2));
    elseif BiDiPhase<0
        minDs(2) = min(BiDiPhase, minDs(2));
    end
    
    ops.yrange = ceil(maxDs(1)):floor(Ly+minDs(1));
    ops.xrange = ceil(maxDs(2)):floor(Lx+minDs(2));
else
    ops.yrange = 1:Ly;
    ops.xrange = 1:Lx;    
end

savepath = sprintf('%s/%s/%s', ops.ResultsSavePath, ops.mouse_name, ops.date);

if ~exist(savepath, 'dir')
    mkdir(savepath)
end
save(sprintf('%s/%s/%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, ...
    ops.mouse_name, ops.date, iplane),  'ops')
%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
