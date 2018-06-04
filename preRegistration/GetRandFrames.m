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
nfmax = max(1, round(ops.NimgFirstRegistration/ntifs)); % number of tiffs to take per file
if nfmax>=2000
    nfmax = 1999;
end

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

% grab random frames from target experiment
if ~isempty(targetImage)
    k = targetImage(1);
    % compute number of frames in tiff
    if abs(nbytes - fs{k}(1).bytes)>1e3
        nbytes = fs{k}(1).bytes;
        nFr = nFrames(fs{k}(1).name);
    end
    startTiff = floor(length(fs{k})/2 - ...
        nplanes*nchannels*ops.NimgFirstRegistration/2/nFr);
    iplane = mod(((startTiff-1)*nFr)/nchannels-1, nplanes) + 1;
    j = startTiff;
    while indx < ops.NimgFirstRegistration
        offset = 0;
        if j == 1
            offset = nchannels*nplanes;
        end
        % compute number of frames in tiff if size different from previous
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFrames(fs{k}(j).name);
        end
        numFr = min(ops.NimgFirstRegistration-indx, nFr/nchannels/nplanes);
        % load red channel if red_align == 1
        if red_align
            ichanset = [offset + nchannels*(nplanes-iplane) + [rchannel;...
                nchannels*nplanes*numFr]; nchannels];
        else
            ichanset = [offset + nchannels*(nplanes-iplane) + [ichannel;...
                nchannels*nplanes*numFr]; nchannels];
        end
        iplane = mod(iplane + nFr/nchannels - 1, nplanes) + 1;
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ichanset(3));
        data = reshape(data, Ly, Lx, nplanes, []);
        
        IMG(:,:,:,indx+(1:size(data,4))) = data;
        indx = indx + size(data,4);
        j = j + 1;
    end
% grab frames from all files
else
    for k = 1:length(ops.SubDirs)
        iplane0 = 1;
        nchannels = ops.nchannels;
        if ~isempty(ops.expts)
            if ismember(ops.expts(k), getOr(ops, 'expred', []))
                nchannels = ops.nchannels_red;
            end
        end
        
        for j = 1:length(fs{k})
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
end
IMG(:,:,:,(1+indx):end) = [];
