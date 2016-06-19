function Fcell = get_signals(ops, iplane)

try
   load(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat',...
       ops.ResultsSavePath,ops.mouse_name, ops.date, iplane, ops.Nk))
catch
   error('Could not find cell detection file \n') 
end

Nk = numel(stat);
Nkpar = ops.Nk;




%% get signals  
[Ly Lx] = size(ops.mimg);

nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F = zeros(Nk, sum(ops.Nframes), 'single');
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
       break; 
    end
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    NT= size(data,3);
    
    data = single(reshape(data, [], NT));
    
    for k = 1:Nk
       ipix = stat(k).ipix; 
       if ~isempty(ipix)
%            F(k,ix + (1:NT)) = stat(k).lambda' * data(ipix,:);       
            F(k,ix + (1:NT)) = mean(data(ipix,:), 1);       
       end
    end    
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);
% F = F(:, 1:ix);

csumNframes = [0 cumsum(ops.Nframes)];
Fcell = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i} = F(:, csumNframes(i) + (1:ops.Nframes(i)));
end

save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, ops.Nk),  'ops', 'res', 'stat', 'stat0', 'res0', 'Fcell', 'clustrules')

