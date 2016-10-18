function db = copy_from_zserver(db, ops)

if ~iscell(db.mouse_name)
    mouse_name{1} = db.mouse_name;
    date{1} = db.date;
    expts{1} = db.expts;
else
    mouse_name = db.mouse_name;
    date = db.date;
    expts = db.expts;
end

for i = 1:length(mouse_name)
    TempDir = fullfile(ops.TempStorage, mouse_name{1}, date{i});
    RootDir = fullfile(ops.RootStorage, mouse_name{i}, date{i});
    if ~exist(TempDir, 'dir');  mkdir(TempDir); end
    for k = 1:length(expts{i})
        SubDirs = num2str(expts{i}(k));
        
        
        fs = dir(sprintf('%s//%s//*.tif', RootDir, SubDirs));
        if ~exist(fullfile(TempDir, SubDirs), 'dir')
            mkdir(fullfile(TempDir, SubDirs));
        end
        
        for j = 1:length(fs)
            if ~exist(fullfile(TempDir, SubDirs, fs(j).name), 'file')
                copyfile(fullfile(RootDir, SubDirs,fs(j).name),...
                    fullfile(TempDir, SubDirs,fs(j).name))
            end
        end
        
%         if ~exist(fullfile(TempDir, SubDirs), 'dir')
%             copyfile(fullfile(RootDir, SubDirs), fullfile(TempDir, SubDirs))
%         end
    end
end

if iscell(db.mouse_name)
   db.mouse_name =  mouse_name{1};
   db.date = date{1};
   db.expts = [];
   for i = 1:length(mouse_name)
      db.expts = [db.expts expts{i}]; 
   end
end