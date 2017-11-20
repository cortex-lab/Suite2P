% non-rigid registration of frames with offsets ds
function [dreg, Valid]= nonrigidRegFrames(data, ops, ds)

orig_class = class(data);

% if ops.useGPU
%     data = gpuArray(single(data));
% end
[Ly, Lx, NT] = size(data);

%%

% % smooth offsets across blocks by xyMask
% nblocks = size(ds,3);
% xyB = zeros(nblocks,2,'single');
% for ib = 1:nblocks
%     xyB(ib,:) = [mean(ops.yBL{ib}) mean(ops.xBL{ib})];
% end
% shiftXY = permute(ds,[3 1 2]);
% 
% dx=interp1(xyB(:,1),shiftXY(:,:,1), [1:Ly]', 'pchip','extrap');
% dy=interp1(xyB(:,1),shiftXY(:,:,2), [1:Ly]', 'pchip','extrap');
% dx = repmat(permute(dx,[1 3 2]),1,Lx,1);
% dy = repmat(permute(dy,[1 3 2]),1,Lx,1);
% dx = round(dx);
% dy = round(dy);
%%
dx = round(ops.xyMask * reshape(squeeze(ds(:,2,:))',size(ops.xyMask,2),[]));
dy = round(ops.xyMask * reshape(squeeze(ds(:,1,:))',size(ops.xyMask,2),[]));

dx = reshape(dx, Ly, Lx, []);
dy = reshape(dy, Ly, Lx, []);

idy = repmat([1:Ly]', 1, Lx);
idx = repmat([1:Lx],  Ly, 1);

dreg = zeros(size(data), orig_class);
Valid = true(Ly, Lx);
for i = 1:NT
    Im = data(:,:,i);    
    
    % apply offsets to indexing
    DX = dx(:,:,i) + idx;
    DY = dy(:,:,i) + idy;
    
    
    % compute valid area of frame
    xyInvalid = DX<0 | DX>Lx-1 | DY<1 | DY>Ly;
    Valid(xyInvalid) = false;
    
    DX(xyInvalid) = 0;
    DY(xyInvalid) = 1;
    
%     DX = mod(DX, Lx);
%     DY = mod(DY-1, Ly) + 1;
%     
    ind = DY + DX * Ly;
    Im = Im(ind);
    dreg(:,:,i) = Im;
    
    
end

