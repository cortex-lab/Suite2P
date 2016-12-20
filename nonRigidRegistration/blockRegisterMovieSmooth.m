function [dreg, Valid, ds]= blockRegisterMovieSmooth(data, ops, ds)

orig_class = class(data);

% if ops.useGPU
%     data = gpuArray(single(data));
% end
[Ly, Lx, NT] = size(data);

%%
numBlocks = ops.numBlocks;
xyMask    = ops.xyMask;

xg = [1:Ly];
xt = [];
for ib = 1:numBlocks
    xt(ib) = round(mean(ops.yBL{ib}));
end

if ops.useGPU
    xg = gpuArray(single(xg));
    xt = gpuArray(single(xt));
end
    
if ops.kriging
    % compute kernels for regression
    % within frame smoothing
    xm  = xt - mean(xt);
    xm2 = sum(xm.^2);
    for j = 1:2
        f   = squeeze(ds(:,j,:));
        fm  = bsxfun(@minus, f, mean(f,2));
        a   = sum(bsxfun(@times, fm, xm),2) / xm2;
        b   = mean(f,2) - a*mean(xt);
        dxg = bsxfun(@times, a, repmat(xg, numel(a), 1));
        dxg = bsxfun(@plus, dxg, b);
        dxg = repmat(dxg',1,1,Lx);
        dxg = permute(dxg,[1 3 2]);
        if j == 1
            dy = dxg;
        else
            dx = dxg;
        end
    end
    clear dxg Kx Kg Kmat;
else
    dx1 = xyMask(:,1:numBlocks) * squeeze(ds(:,2,:))';
    dy1 = xyMask(:,1:numBlocks) * squeeze(ds(:,1,:))';
    
    dx1 = reshape(dx1, Ly, Lx, []);
    dy1 = reshape(dy1, Ly, Lx, []);
    dx = dx1;
    dy = dy1;
end
ds = [dy(:,1,:) dx(:,1,:)];

%keyboard;
% across frame smoothing (only good in single plane imaging)
if 0 
    % temporal smoothing in y (over line scanning)
    tsm = ops.regSmooth;
    % first lines smoothed with last lines on prev frame
    dx(ops.yBL{1}, 2:NT) = (1-tsm) * dx(ops.yBL{1}, 2:NT) + tsm * dx(ops.yBL{numBlocks}, 1:NT-1);
    % last lines smoothed with first lines on next frame
    dx(ops.yBL{numBlocks}, 1:NT-1) = (1-tsm) * dx(ops.yBL{numBlocks}, 1:NT-1) + tsm * dx(ops.yBL{1}, 2:NT);
end

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

