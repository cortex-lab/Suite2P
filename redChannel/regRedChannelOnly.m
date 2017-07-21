function mimgR = regRedChannelOnly(ops)

% numPlanes = length(ops.planesToProcess);

%% build file list with red channel

if (isfield(ops, 'SubDirsRed') && ~isempty(ops.SubDirsRed))
   subDirsRed = ops.SubDirsRed;
else
    if (isfield(ops, 'expred') && ~isempty(ops.expred))
        for i = 1:length(ops.expred)
            subDirsRed{i} = sprintf('%d', ops.expred(i));
        end
    else
        warning('could not find red channel info, returning...')
        return;
    end 
end
%%
% build file list
fsrootRED = [];

for j = 1:length(subDirsRed)
    fsrootRED{j} = dir(fullfile(ops.RootDir, subDirsRed{j}, '*.tif'));
    
    for k = 1:length(fsrootRED{j})
        fsrootRED{j}(k).name = fullfile(ops.RootDir, subDirsRed{j}, fsrootRED{j}(k).name);         
    end
end

%%
fsRED = [];
for j = 1:length(fsrootRED)
    fsRED = [fsRED {fsrootRED{j}.name}];
end


%%
D = [];
for k = 1:length(fsRED)    
    data = loadFramesBuff(fsRED{k});
    
    data = cat(3, D, data);    
end


%%
clear mimgR
for j = 1:ops.nplanes
    i0 = ops.nchannels_red + (j-1)*ops.nchannels_red; 
    %mimgR(:,:,j) = mean(data(:,:,i0:ops.nplanes*ops.nchannels:end), 3);
    datj = data(:,:,i0:ops.nplanes*ops.nchannels_red:end);
    ops1  = AlignIterativeKriging(datj, ops);
    mimgR(:,:,j) = ops1.mimg;    
end

%%