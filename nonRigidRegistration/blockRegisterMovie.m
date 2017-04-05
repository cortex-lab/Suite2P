% non-rigid registration of frames with offsets ds
function [dreg, Valid]= blockRegisterMovie(data, xyMask, ds)

orig_class = class(data);

% if ops.useGPU
%     data = gpuArray(single(data));
% end
[Ly, Lx, NT] = size(data);

%%

% smooth offsets across blocks by xyMask
dx = round(xyMask * squeeze(ds(:,2,:))');
dy = round(xyMask * squeeze(ds(:,1,:))');

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

