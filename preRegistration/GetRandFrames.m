% gets random frames (in total NimgFirstRegistration)
% frames are used for initializing registration
function IMG = GetRandFrames(fs, ops)

nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

% if not empty, generate target image from specified experiment (useful if there is slow drift in data)
targetImage        = getOr(ops, {'targetImage'}, []); 

ntifs = sum(cellfun(@(x) numel(x), fs)); % total number of tiff files

% size of images
nbytes = fs{1}(1).bytes;
Info0 = imfinfo(fs{1}(1).name);
nFr = length(Info0);
Ly = Info0(1).Height;
Lx = Info0(1).Width;
ops.Lx = Lx;
ops.Ly = Ly;

indx = 0;
IMG = zeros(Ly, Lx, nplanes, ops.NimgFirstRegistration, 'single');

if ~isempty(targetImage) % grab frames around half time during target experiment
    expInds = targetImage(1);
    startTiff = floor(length(fs{expInds})/2 - ...
        nplanes*nchannels*ops.NimgFirstRegistration/2/nFr);
    tiffInds = startTiff:length(fs{expInds});
    iplaneK = 1 - (startTiff-1)*nFr; % to calculate index of first frame belonging to plane 1 in startTiff
    nfmax = floor(nFr/nplanes/nchannels)-1;
else
    expInds = 1:length(ops.SubDirs);
    tiffInds = 1:length(fs{expInds});
    iplaneK = 1;
    nfmax = max(1, round(ops.NimgFirstRegistration/ntifs)); % number of tiffs to take per file
    if nfmax>=2000
        nfmax = 1999;
    end
end

for k = expInds
    iplane0 = iplaneK;
    nchannels = getOr(ops, {'nchannels'}, 1);
    if ~isempty(ops.expts)
        if ismember(ops.expts(k), getOr(ops, 'expred', []))
            nchannels = ops.nchannels_red;
        end
    end
    
    for j = tiffInds
        % compute number of frames in tiff if size different from previous
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFrames(fs{k}(j).name);
        end
        if j==1
            offset = nchannels*nplanes;
        else
            offset = 0;
        end
        % check if there are enough frames in the tiff
        if nFr<(nchannels*nplanes*nfmax+offset)
            continue;
        end
        
        iplane0 = mod(iplane0-1, nplanes) + 1;
        % load red channel if red_align == 1
        if red_align
            ichanset = [offset + nchannels*(iplane0-1) + [rchannel;...
                nchannels*nplanes*nfmax]; nchannels];
        else
            ichanset = [offset + nchannels*(iplane0-1) + [ichannel;...
                nchannels*nplanes*nfmax]; nchannels];
        end
        iplane0 = iplane0 - nFr/nchannels;
        data = loadFramesBuff(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
        data = reshape(data, Ly, Lx, nplanes, []);
        
        IMG(:,:,:,indx+(1:size(data,4))) = data;
        indx = indx + size(data,4);
        
        if indx>ops.NimgFirstRegistration
            break;
        end
    end
    
    if indx>ops.NimgFirstRegistration
        break;
    end
end

IMG(:,:,:,(1+indx):end) = [];
