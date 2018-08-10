function write_reg_to_tiff(fid, ops, iplane, isRED, sm)

Ly = ops.Ly;
Lx = ops.Lx;
bitspersamp = 16;

frewind(fid);
numPartitions = sm.getNumPartitions();
for k = 1:numPartitions
    sm.seekToPartition(k);
    
    ix = 0;
    nframesleft = ops.Nframes(k);

    while nframesleft > 0
        ix = ix + 1;
        nfrtoread = min(nframesleft, 2000);
        data = fread(fid,  Ly*Lx*nfrtoread, 'int16');
        nframesleft = nframesleft - nfrtoread;
        data = reshape(data, Ly, Lx, []);

        if isRED
            fname = sm.registrationTiffPath(iplane, ix, 'red');
        else
            fname = sm.registrationTiffPath(iplane, ix);
        end

        TiffWriter(int16(data), fname, bitspersamp);
    end
end
