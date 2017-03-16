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

root = ops.ResultsSavePath;
fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
if exist(fullfile(root, fregops), 'file')
    load(fullfile(root, fregops))
else
    ops1 = cell(ops.nplanes, 1);
    for j = 1:ops.nplanes
        
        fname = sprintf('regops_%s_%s_plane%d.mat', ops.mouse_name, ops.date, j);
        dat = load(fullfile(root, fname));
        ops1{j} = dat.ops;
        ops1{j}.useGPU = ops.useGPU;
    end
end

%
ntf0 = 0;
numPlanes = ops.nplanes;
%iplane0 = 1:1:ops.nplanes;

totFrames=0;
for k = 1:length(fsRED)
    %iplane0 = mod(iplane0-1, numPlanes) + 1;
    startPlane=((mod(totFrames+1,ops.nplanes*2)-1)/2)+1;
    
    nFr = nFramesTiff(fsRED(k).name);
    totFrames=totFrames+nFr;
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
    
    
    for iPlane=1:ops.nplanes
      
        idx0=mod((ops.nplanes-startPlane+iPlane)*2,ops.nplanes*2);
        planesG=(idx0+1):(2*ops.nplanes):nFr;
        planesR=(idx0+2):(2*ops.nplanes):nFr;
        dataG0=data(:,:,planesG);
        dataR0=data(:,:,planesR);
        
        BiDiPhase = ops1{iPlane}.BiDiPhase;
        if abs(BiDiPhase) > 0
            yrange = 2:2:Ly;
            if BiDiPhase>0
                dataG0(yrange,(1+BiDiPhase):Lx,:,:) = dataG0(yrange, 1:(Lx-BiDiPhase),:,:);
                dataR0(yrange,(1+BiDiPhase):Lx,:,:) = dataR0(yrange, 1:(Lx-BiDiPhase),:,:);
            else
                dataG0(yrange,1:Lx+BiDiPhase,:,:)   = dataG0(yrange, 1-BiDiPhase:Lx,:,:);
                dataR0(yrange,1:Lx+BiDiPhase,:,:)   = dataR0(yrange, 1-BiDiPhase:Lx,:,:);
            end
        end
        
        [ds, ~]  = registration_offsets(dataG0, ops1{iPlane}, 0);
        if k==1
            ds(1,:) = 0;
        end
        dataR       = ...
            register_movie(dataR0, ops1{iPlane}, ds);
        dataG     = ...
            register_movie(dataG0, ops1{iPlane}, ds);
        
        
        mimgR(:,:,iPlane) = mimgR(:,:,iPlane) + mean(dataR, 3);
        mimgG(:,:,iPlane) = mimgG(:,:,iPlane) + mean(dataG, 3);
        
    end
    
    ntf0 = ntf0 + 1;
    
    %iplane0 = iplane0 - nFr/ops.nchannels_red;
    fprintf('processing tiff %d/%d\n',k,length(fsRED))
end

mimgR = mimgR/ntf0;
mimgG = mimgG/ntf0;