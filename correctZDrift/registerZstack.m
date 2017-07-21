% registers Z-stack using kriging phase correlation
% saves Z-stack to ops.ZstackSavePath
function registerZstack(ops)

zPath = fullfile(ops.RootDir, num2str(ops.Zexpt));
fs = dir(fullfile(zPath, '*.tif'));
for i = 1:length(fs)
   fs(i).name = fullfile(zPath, fs(i).name);
end

nFr0 = nFrames(fs(1).name);
%
figure;
for j = 1:ops.Zplanes
    ifirst = 0; % which channel we want to image (0 for first, 1 for second)
    nFr = nFr0;
    I = [];
    for i = 1:length(fs)
        if i==length(fs)
            nFr = nFrames(fs(i).name);
        end
        i0 = mod(ifirst+ops.Zchannels*j-1, ops.Zchannels*ops.Zplanes) + 1;
        I = cat(3, I, loadFramesBuff(fs(i).name, i0, nFr, ops.Zplanes * ops.Zchannels));
        
        ifirst = ifirst - nFr;
        
        if j==1
            Mimg = zeros(size(I,1), size(I,2), ops.Zplanes*ops.Zchannels);
        end
    end
    ops.Ly = size(I,1);
    ops.Lx = size(I,2);
    ops = AlignIterativeKriging(single(I), ops);

    Mimg(:,:,j) = ops.mimg;
    
    clf
    imagesc(Mimg(:,:,j));
    title(sprintf('Z-stack: plane %d', j));
    drawnow;
end

%
Mimg0 = single(Mimg);
%
fname = sprintf('stack_%s_%s.mat', ops.mouse_name, ops.date);
save(fullfile(ops.ZstackSavePath, fname), 'Mimg0')

