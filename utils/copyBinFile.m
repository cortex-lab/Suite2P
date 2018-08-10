function copyBinFile(ops1, sm)
    for i = 1:size(ops1, 1)
        for j = 1:size(ops1, 2)
            fid = fopen(ops1{i, j}.RegFile, 'r');

            fileName = sm.registrationFileOnSlowStorage(i);
            fidCopy = fopen(fileName, 'w');

            sz = ops1{i,j}.Lx * ops1{i,j}.Ly;
            parts = ceil(sum(ops1{i, j}.Nframes) / 2000);
            for p = 1:parts
                toRead = 2000;
                if p == parts
                    toRead = sum(ops1{i,j}.Nframes) - 2000 * (parts-1);
                end
                data = fread(fid,  sz*toRead, '*int16');
                fwrite(fidCopy, data, class(data));
            end
            fclose(fidCopy);
            fclose(fid);
        end
    end
end