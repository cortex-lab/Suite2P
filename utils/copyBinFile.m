function copyBinFile(ops1)

fid    = cell(size(ops1));
for i = 1:size(ops1,1)
    for j = 1:size(ops1,2)
        fid{i,j}           = fopen(ops1{i,j}.RegFile, 'r');
        folder = fullfile(ops1{i,j}.RegFileBinLocation, ops1{i,j}.mouse_name, ...
            ops1{i,j}.date, ops1{i,j}.CharSubDirs);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        fidCopy = fopen(fullfile(folder, ...
            sprintf('%s_%s_%s_plane%d.bin', ops1{i,j}.mouse_name, ops1{i,j}.date, ...
            ops1{i,j}.CharSubDirs, i)), 'w');
        sz = ops1{i,j}.Lx * ops1{i,j}.Ly;
        parts = ceil(sum(ops1{i,j}.Nframes) / 2000);
        for p = 1:parts
            toRead = 2000;
            if p == parts
                toRead = sum(ops1{i,j}.Nframes) - 2000 * (parts-1);
            end
            data = fread(fid{i,j},  sz*toRead, '*int16');
            fwrite(fidCopy, data, class(data));
        end
        fclose(fidCopy);
        fclose(fid{i,j});
    end
end