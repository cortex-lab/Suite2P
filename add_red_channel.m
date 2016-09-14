function add_red_channel(ops, mimgR,mimgG)


for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    try
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        while isfield('dat', 'dat')
            dat = dat.dat;
        end
        dat.ops.mimgRED = mimgR(:,:,iplane);
        if ~isempty(mimgG)
            dat.ops.mimgGREEN = mimgG(:,:,iplane);
        end
        save(fname, 'dat')
    catch
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        dat.ops.mimgRED = mimgR(:,:,iplane);
        if ~isempty(mimgG)
            dat.ops.mimgGREEN = mimgG(:,:,iplane);
        end
        save(fname, '-struct', 'dat')
    end
    
    fname = sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
    dd = load(fname);
    dd.ops.mimgRED = mimgR(:,:,iplane);
    if ~isempty(mimgG)
        dd.ops.mimgGREEN = mimgG(:,:,iplane);
    end
    save(fname, '-struct', 'dd')
end
