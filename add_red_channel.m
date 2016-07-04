function add_red_channel(ops, mimgR)


for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    try
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        while isfield('dat', 'dat')
            dat = dat.dat;
        end
        dat.ops.mimgRED = mimgR(:,:,iplane);
        save(fname, 'dat', '-v7.3')
    catch
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        dat.ops.mimgRED = mimgR(:,:,iplane);
        save(fname, '-struct', 'dat', '-v7.3')
    end
    
    
    fname = sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
    dd = load(fname);
    dd.ops.mimgRED = mimgR(:,:,iplane);
    save(fname, '-struct', 'dd')
end
