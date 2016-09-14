function ops = getPVimage(ops)

root = fullfile(ops.RootDir, num2str(ops.expred));
fs = dir(fullfile(root, '*.tif'));

nmaximg = 500;

I1 = [];
I2 = [];
for j = 1:numel(fs)
    mov = img.loadFrames(fullfile(root, fs(j).name), 1,2*nmaximg,1);
    I2 = cat(3, I2, mov(:,:,2:2:end));
    I1 = cat(3, I1, mov(:,:,1:2:end));
    if size(I2,3)>nmaximg
       break; 
    end
end
I2 = I2(:,:,1:nmaximg);
I1 = I1(:,:,1:nmaximg);
I1 = single(I1);
I2 = single(I2);

ops1    = align_iterative(I2, ops);
I1reg   = register_movie(I1, ops1, ops1.dsprealign);
mI1     = mean(I1reg,3);

[dsnew, Corr]  = registration_offsets(mI1, ops, 0);
PVreg   = single(register_movie(ops1.mimg, ops1, dsnew));


ops.imgPV = PVreg;