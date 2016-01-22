function Fcell = get_signals_and_neuropil(ops, iplane)

try
    load(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat',...
        ops.ResultsSavePath,ops.mouse_name, ops.date, iplane, ops.Nk))
catch
    error('Could not find cell detection file \n')
end
%%

Nk        = numel(stat);   % 1277, Nk > Nkparents, all candidate ROIs are > NkParents
NkParents = ops.Nk; % 650
[LyU, LxU] = size(ops.mimg);
LyR=length(ops.yrange);
LxR=length(ops.xrange);

iscellID=(NkParents+1:Nk);
nCells1=length(iscellID);


allField=zeros(LyR,LxR);

cellFields=nan(nCells1,LyR,LxR);
for jCell=1:nCells1
    iCell=iscellID(jCell);
    ipix=stat(iCell).ipix;
    temp=zeros(LyR,LxR);
    temp(ipix)=ones(1,length(ipix));
    cellFields(jCell,:,:)=temp;
    allField=allField+jCell.*temp;
end

opt.inNeurop=3; %fixed inner diameter of the neuropil mask donut
opt.outNeurop=45; %radius of Neuropil fixed at 45um
opt.zoomMicro=2; %fixed zoomMicro
opt.microID='b';
opt.totPixels=512; %fixed number of pixel


um2pix=infoPixUm(opt.totPixels,opt.zoomMicro,opt.microID);
xPU=um2pix.xPU;
yPU=um2pix.yPU;
neuropMasks=createNeuropilMasks(cellFields,allField,xPU,yPU,opt);

%Note: this part is quite slow
fprintf('Slow part\n')
tic
mCell=0;
for k=iscellID
    mCell=mCell+1;
    temp=squeeze(neuropMasks(mCell,:,:));
    stat(k).ipix_neuropil=find(temp==1);
end
toc


% his function needs to add a field of pixel indices to stat called ipix_neuropil, only for ROIs from NkParents+1 to Nk
%stat = add_neuropil_ROI(stat);

%% get signals


nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F = zeros(Nk, sum(ops.Nframes), 'single');
Fneu = zeros(Nk-NkParents, sum(ops.Nframes), 'single');

while 1
    data = fread(fid,  LyU*LxU*nimgbatch, '*int16');
    if isempty(data)
        break;
    end
    data = reshape(data, LyU, LxU, []);
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
        
        if k>NkParents
            ipix_neuropil= stat(k).ipix_neuropil;
            if ~isempty(ipix_neuropil)
                Fneu(k-NkParents,ix + (1:NT)) = mean(data(ipix_neuropil,:), 1);
            end
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
FcellNeu = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i} 	= F(:, csumNframes(i) + (1:ops.Nframes(i)));
    FcellNeu{i} = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end

save(sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date, iplane, ops.Nk),  'ops', 'res', 'stat', 'stat0', 'res0', 'Fcell', 'FcellNeu')

