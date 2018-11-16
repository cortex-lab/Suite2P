% gets random frames (in total NimgFirstRegistration)
% frames are used for initializing registration
function IMG = GetFirstFrames(fs, ops)

nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

% if not empty, generate target image from specified experiment (useful if there is slow drift in data)
targetImage        = getOr(ops, {'targetImage'}, []); 


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

for k = 1:length(ops.SubDirs)
    iplane0 = 1;
    if ismember(ops.expts(k), getOr(ops, 'expred', []))
        nchannels = ops.nchannels_red;
    else
        nchannels = ops.nchannels;
    end
    for j = 1:length(fs{k})
        % compute number of frames in tiff if size different from previous
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFrames(fs{k}(j).name);
        end
        nrez = ops.NimgFirstRegistration - indx;
        if j==1
            offset = nchannels*nplanes;
            nFrMax = min(nFr-nplanes * nchannels, nrez* nplanes * nchannels);
        else
            offset = 0;
            nFrMax = min(nFr, nrez * nplanes * nchannels);
        end
                        
        iplane0 = mod(iplane0-1, nplanes) + 1;
        % load red channel if red_align == 1
                
        if red_align
            ichanset = [offset + nchannels*(iplane0-1) + [rchannel;...
                nFrMax]; nchannels];
        else
            ichanset = [offset + nchannels*(iplane0-1) + [ichannel;...
                nFrMax]; nchannels];
        end
        iplane0 = iplane0 - nFr/nchannels;
        
        data = loadFramesBuff(fs{k}(j).name, ichanset(1),ichanset(2), ichanset(3));
        data = reshape(data, Ly, Lx, nplanes, []);
        
        IMG(:,:,:,indx+(1:size(data,4))) = data;
        indx = indx + size(data,4);
        
        if indx>=ops.NimgFirstRegistration
            break;
        end
    end
end

IMG(:,:,:,(1+indx):end) = [];
