% function to load and view the registered binaries output from Suite2P

% db file
make_db_MP024;
% which experiment to view
iexp = 5;

db = db(iexp);
ops.ResultsSavePath = 'D:/DATA/F/';
CharSubDirs = '';
for k = 1:length(db.expts)
    CharSubDirs    = [CharSubDirs num2str(db.expts(k)) '_'];
end
CharSubDirs = CharSubDirs(1:end-1);

ops.ResultsSavePath = sprintf('%s//%s//%s//%s//', ops.ResultsSavePath, db.mouse_name, db.date, ...
    CharSubDirs);
regpath = sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, db.mouse_name, db.date);
load(regpath);

%%
iplane = 5; % which plane of the recording
ops = ops1{iplane};
fid  = fopen(ops.RegFile, 'r'); % opens the registered binary
[Ly Lx] = size(ops.mimg1);    % size of binary in x,y

clf;
nt0 = 0;
nimgbatch = 2000; % how many frames you can load at a time (will depend on RAM)
while 1
  data = fread(fid, Ly*Lx*nimgbatch, '*int16');
  if isempty(data)
    break;
  end
  data = reshape(data, Ly, Lx, []);
  data = data(ops.yrange, ops.xrange, :);
  NT   = size(data,3);
  for j = 1:NT
    imagesc(data(:,:,j),[0 2000]);%,[1000 3000]);
    title(sprintf('frame %d', j+nt0));
    axis square;
    drawnow;
    
    pause(0.25);
  end
  nt0 = nt0 + NT;
end
fclose all;
