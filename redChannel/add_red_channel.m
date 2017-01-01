function add_red_channel(ops, mimgR,mimgG)


root = ops.ResultsSavePath;
fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
flag = 0;
if exist(fullfile(root, fregops), 'file')
    dd = load(fullfile(root, fregops));
    flag = 1;
end

for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    try
        try
            fname = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
            dat = load(fname);
            while isfield(dat, 'dat')
                dat = dat.dat;
            end
            dat.ops.mimgRED = mimgR(:,:,iplane);
            if ~isempty(mimgG)
                dat.ops.mimgGREEN = mimgG(:,:,iplane);
            end
            save(fname, 'dat')
        catch
            fname = sprintf('%s/F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
            dat = load(fname);
            while isfield(dat, 'dat')
                dat = dat.dat;
            end
            dat.ops.mimgRED = mimgR(:,:,iplane);
            if ~isempty(mimgG)
                dat.ops.mimgGREEN = mimgG(:,:,iplane);
            end
            save(fname, 'dat')
        end
    catch
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        dat.ops.mimgRED = mimgR(:,:,iplane);
        if ~isempty(mimgG)
            dat.ops.mimgGREEN = mimgG(:,:,iplane);
        end
        save(fname, '-struct', 'dat')
    end
    
    if flag 
        dd.ops1{i}.mimgRED = mimgR(:,:,iplane);
        if ~isempty(mimgG)
            dd.ops1{i}.mimgGREEN = mimgG(:,:,iplane);
        end
    else
        fname = sprintf('%s/regops_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
        dd = load(fname);
        dd.ops.mimgRED = mimgR(:,:,iplane);
        if ~isempty(mimgG)
            dd.ops.mimgGREEN = mimgG(:,:,iplane);
        end
        save(fname, '-struct', 'dd')
    end
end

if flag
    save(fullfile(root, fregops), '-struct', 'dd')
end
