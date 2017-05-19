function dat = get_signals_and_neuropil(opt, iplane, doSave)

if ~isfield(opt,'zoomMicro')
    opt.zoomMicro=2; %fixed zoomMicro
end
if ~isfield(opt,'inNeurop')
    opt.inNeurop=3; %fixed inner diameter of the neuropil mask donut
end
if ~isfield(opt,'outNeurop')
    opt.outNeurop=45; %radius of Neuropil fixed at 45um
end
if ~isfield(opt,'microID')
    opt.microID='b'; %microscope identity
end
if ~isfield(opt, 'processed') % use processed data in F_...._proc (generated in gui2P)
    opt.processed = 0;
end
if ~isfield(opt, 'useSVD') % redo calculation of signal and neuropil based on
                         % saved SVD components instead of temporary
                         % bin-file
    opt.useSVD = 0;
end
if ~isfield(opt, 'getSignal') % calculate and save neural signal
    opt.getSignal = 1;
end
if ~isfield(opt, 'getNeuropil') % calculate and save neuropil
    opt.getNeuropil = 1;
end
if ~isfield(opt, 'newFile') % save new file '<name>_new.mat', otherwise
                            % existing file is overwritten
    opt.newFile = 0;
end
if nargin < 3
    doSave = 1;
end

filenames = dir(sprintf('%s/F_%s_%s_plane%d*.mat',...
    opt.ResultsSavePath, opt.mouse_name, opt.date, iplane));
filenames = {filenames.name};
if isfield(opt, 'Nk')
    ind = cellfun(@strfind, filenames, ...
        repmat({['_Nk' num2str(opt.Nk)]}, size(filenames)), ...
        'UniformOutput', false);
    ind = ~cellfun(@isempty, ind);
    filenames = filenames(ind);
end
ind = cellfun(@strfind, filenames, repmat({'_proc'}, size(filenames)), ...
    'UniformOutput', false);
ind = ~cellfun(@isempty, ind);
if opt.processed == 1
    filenames = filenames(ind);
else
    filenames = filenames(~ind);
end
if length(filenames) > 1
    fprintf(['WARNING: several files for dataset exist\n' ...
        '%s <- using\n'], filenames{1})
    fprintf('%s\n', filenames{2:end})
elseif isempty(filenames)
    error('Could not find cell detection file \n')
end

data = load(fullfile(opt.ResultsSavePath, filenames{1}));
if opt.processed == 1
    data = data.dat;
end
%%
ops = data.ops;
Nk       = numel(data.stat); % all ROIs

useCells = 1:Nk;
if opt.processed == 1
    useCells = find([data.stat.iscell]);
end
[LyU, LxU] = size(ops.mimg);
LyR=length(ops.yrange);
LxR=length(ops.xrange);

allField=zeros(LyR,LxR);

cellFields=nan(length(useCells),LyR,LxR);
exclude = [];
for jCell=1:length(useCells)
    iCell=useCells(jCell);
    % use all pixels of ROI to include into allFields (used to define
    % neuropil mask)
    ipix=data.stat(iCell).ipix;
    temp=zeros(LyR,LxR);
    temp(ipix)=ones(1,length(ipix));
    allField=allField+jCell.*temp;
    if isfield(data.stat, 'isoverlap')
        % if overlap of this ROI with other ROIs is larger than 2/3s of its
        % size, exclude this ROI
        overlap = data.stat(iCell).isoverlap;
        if sum(overlap) > .66 * length(ipix)
            exclude(end+1) = jCell;
            data.stat(iCell).iscell = 0;
            continue
        end
        % exclude pixels from ROI that overlap with other ROIs to determine
        % signal of this ROI
        ipix(overlap) = [];
    end
    temp=zeros(LyR,LxR);
    temp(ipix)=ones(1,length(ipix));
    cellFields(jCell,:,:)=temp;
end
cellFields(exclude,:,:) = [];
useCells(exclude) = [];
opt.totPixels=LxU;

um2pix=infoPixUm(opt.totPixels,opt.zoomMicro,opt.microID);
xPU=um2pix.xPU;
yPU=um2pix.yPU;
if opt.getNeuropil
    neuropMasks=createNeuropilMasks(cellFields,allField,xPU,yPU,opt);

    mCell=0;
    for k=useCells
        mCell=mCell+1;
        temp=squeeze(neuropMasks(mCell,:,:));
        data.stat(k).ipix_neuropil=find(temp==1);
    end
end

%% get signals and neuropil
if opt.useSVD == 0
    nimgbatch = 2000;
    ix = 0;
    fclose all;
    fid = fopen(ops.RegFile, 'r');
    
    tic
    F = NaN(Nk, sum(ops.Nframes), 'single');
    Fneu = NaN(Nk, sum(ops.Nframes), 'single');
    while 1
        mov = fread(fid,  LyU*LxU*nimgbatch, '*int16');
        if isempty(mov)
            break;
        end
        mov = reshape(mov, LyU, LxU, []);
        mov = mov(ops.yrange, ops.xrange, :);
        mov = single(mov);
        NT= size(mov,3);
        
        mov = single(reshape(mov, [], NT));
        
        for k = 1:length(useCells)
            iCell = useCells(k);
            if opt.getSignal
                ipix = data.stat(iCell).ipix;
                if isfield(data.stat, 'isoverlap')
                    ipix(data.stat(iCell).isoverlap) = [];
                end
                if ~isempty(ipix)
                    % F(k,ix + (1:NT)) = stat(k).lambda' * data(ipix,:);
                    F(iCell,ix + (1:NT)) = mean(mov(ipix,:), 1);
                end
            end
            if opt.getNeuropil
                ipix_neuropil= data.stat(iCell).ipix_neuropil;
                if ~isempty(ipix_neuropil)
                    Fneu(iCell,ix + (1:NT)) = mean(mov(ipix_neuropil,:), 1);
                end
            end
        end
        
        ix = ix + NT;
        if rem(ix, 3*NT)==0
            fprintf('Frame %d done in time %2.2f \n', ix, toc)
        end
    end
    fclose(fid);
    % F = F(:, 1:ix);
    csumNframes = [0 cumsum(ops.Nframes)];
    Fcell = cell(1, length(ops.Nframes));
    FcellNeu = cell(1, length(ops.Nframes));
    for i = 1:length(ops.Nframes)
        Fcell{i} 	= F(:, csumNframes(i) + (1:ops.Nframes(i)));
        FcellNeu{i} = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
    end
else % opt.useSVD == 1
    ind = strfind(filenames{1}, '_Nk');
    fname = ['SVD' filenames{1}(2:ind-1) '.mat'];
    svd = load(fullfile(ops.ResultsSavePath, fname), 'U', 'Vcell');
    U1 = reshape(svd.U, [], size(svd.U,3));
    Ucell = zeros(length(useCells), size(U1,2));
    UcellNeu = zeros(length(useCells), size(U1,2));
    Fcell = cell(size(svd.Vcell));
    FcellNeu = cell(size(svd.Vcell));
    for iCell = 1:length(useCells)
        if opt.getSignal
            ipix = data.stat(useCells(iCell)).ipix;
            if ~isempty(ipix)
                Ucell(iCell,:) = mean(U1(ipix, :), 1);
            end
        end
        if opt.getNeuropil
            ipix_neuropil= data.stat(useCells(iCell)).ipix_neuropil;
            if ~isempty(ipix_neuropil)
                UcellNeu(iCell,:) = mean(U1(ipix_neuropil,:), 1);
            end
        end
    end
    for iExp = 1:length(svd.Vcell)
        if opt.getSignal
            F = NaN(Nk, size(svd.Vcell{iExp}, 2), 'single');
            F(useCells,:) = Ucell * svd.Vcell{iExp};
            Fcell{iExp} = F;
        end
        if opt.getNeuropil
            F = NaN(Nk, size(svd.Vcell{iExp}, 2), 'single');
            F(useCells,:) = UcellNeu * svd.Vcell{iExp};
            FcellNeu{iExp} = F;
        end
    end
end
if opt.processed
    if opt.getSignal
        data.Fcell = Fcell;
    end
    if opt.getNeuropil
        data.FcellNeu = FcellNeu;
    end
else
    if opt.getSignal
        data.Fcell = Fcell;
    end
    if opt.getNeuropil
        data.FcellNeu = FcellNeu;
    end
end

dat = data;
dat.opsNpil = opt;
if doSave == 0
    return
end
if opt.processed == 1
    if opt.newFile == 1
        save(fullfile(opt.ResultsSavePath, [filenames{1}(1:end-4) ...
            '_new.mat']), 'dat')
    else
        save(fullfile(opt.ResultsSavePath, filenames{1}), 'dat');
    end
else
    if opt.newFile == 1
        save(fullfile(opt.ResultsSavePath, [filenames{1}(1:end-4) ...
            '_new.mat']), '-struct', 'dat');
    else
        save(fullfile(opt.ResultsSavePath, filenames{1}), '-struct', 'dat');
    end
end
% save(filenames{1},  'ops', 'res', 'stat', 'stat0', 'res0', 'Fcell', 'FcellNeu')

