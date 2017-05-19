function add_red_channel_sourcery(ops)

load(fullfile(ops.ResultsSavePath,'redchannel.mat'));

root = ops.ResultsSavePath;
fregops =  sprintf('regops_%s_%s.mat', ops.mouse_name, ops.date);
flag = 0;
if exist(fullfile(root, fregops), 'file')
    dd = load(fullfile(root, fregops));
    flag = 1;
end

fclose all;
for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    
    fname = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
    dat = load(fname);
    while isfield(dat, 'dat')
        dat = dat.dat;
    end
    dat.ops.mimgRED = mimgR(:,:,iplane);
    if ~isempty(mimgG)
        dat.ops.mimgGREEN = mimgG(:,:,iplane);
    end
    
    varinfo = whos('dat');
    if varinfo.bytes >= 2^31
        save(fname, '-v7.3', '-struct', 'dat');
    else
        save(fname, '-struct', 'dat');
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
        save(fname, '-v7.3', '-struct', 'dd')
    end
end

if flag
    save(fullfile(root, fregops), '-struct', 'dd')
end
