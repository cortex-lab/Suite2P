function IMG = GetRandFrames(fs, ops)

nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

ntifs = sum(cellfun(@(x) numel(x), fs));
nfmax = max(1, round(ops.NimgFirstRegistration/ntifs));
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

for k = 1:length(ops.SubDirs)
    iplane0 = 1;
    for j = 1:length(fs{k})
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFrames(fs{k}(j).name);
        end
        if j==1
            offset = nchannels*nplanes;
        else
            offset = 0;
        end
        
        if nFr<(nchannels*nplanes*nfmax+offset)
            continue;
        end
        
        iplane0 = mod(iplane0-1, nplanes) + 1;
        
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
IMG =  IMG(:,:,:,1:indx);
