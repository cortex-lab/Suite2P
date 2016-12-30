function [dreg, Valid, ds]= blockRegisterMovieSmooth(data, ops, ds)

ops.regSmooth  = getOr(ops, {'regSmooth'}, 0);
ops.quadBlocks = getOr(ops, {'quadBlocks'}, 1);


orig_class = class(data);

% if ops.useGPU
%     data = gpuArray(single(data));
% end
[Ly, Lx, NT] = size(data);

%%
numBlocks = ops.numBlocks;
xyMask    = ops.xyMask;

xg = [1:Ly]';
xt = [];
for ib = 1:numBlocks
    xt(ib) = round(mean(ops.yBL{ib}));
end
xt = xt';

if ops.useGPU
    xt = gpuArray(single(xt));
end
    
if ops.quadBlocks    
    % within frame smoothing
    xm  = xt - mean(xt);
    xm2 = sum(xm.^2);
    for j = 1:2
        f   = squeeze(ds(:,j,:))';
        dxg = fitQuad(xt, f, xg);
        dxg = repmat(dxg,1,1,Lx);
        dxg = permute(dxg,[1 3 2]);
        if j == 1
            dy = dxg;
        else
            dx = dxg;
        end
    end
    clear dxg xt;
else
    dx1 = xyMask(:,1:numBlocks) * squeeze(ds(:,2,:))';
    dy1 = xyMask(:,1:numBlocks) * squeeze(ds(:,1,:))';
    
    dx1 = reshape(dx1, Ly, Lx, []);
    dy1 = reshape(dy1, Ly, Lx, []);
    dx = dx1;
    dy = dy1;
end
ds = [dy(:,1,:) dx(:,1,:)];

dx = round(dx);
dy = round(dy);



idy = repmat([1:Ly]', 1, Lx);
idx = repmat([1:Lx],  Ly, 1);

dreg = zeros(size(data), orig_class);
Valid = true(Ly, Lx);
for i = 1:NT
    Im = data(:,:,i);    
    
    DX = dx(:,:,i) + idx;
    DY = dy(:,:,i) + idy;
    
    
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

