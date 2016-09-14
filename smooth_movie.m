function dat = smooth_movie(dat, ops)


smooth = reshape(ops.smooth_time_space, 1, []);
switch length(smooth)
    case 1
        smDims = 3;
        smSigs = smooth;
        if smooth <= 0
            smooth = [];
        end
    case 2
        smSigs = smooth([2 2 1]);
        smDims = find(smSigs > 0);
        smSigs = smSigs(smDims);
        if isempty(smDims)
            smooth = [];
        end
    case 3
        smSigs = smooth([2 3 1]);
        smDims = find(smSigs > 0);
        smSigs = smSigs(smDims);
        if isempty(smDims)
            smooth = [];
        end
end
if ~isempty(smooth)
    dat = my_conv2(double(dat), smSigs, smDims);
    dat = int16(dat);
end