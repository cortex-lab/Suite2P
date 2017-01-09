% this function registers data frames by shifts ds and
% returns the registered frames as dreg and the valid
% x-y limits as Valid

function [dreg, Valid, ds0]= blockRegisterMovieSmooth(data, ops, ds)

% whether or not to fit quadratic to ds across a frame
% for interpolation from numBlocks -> Ly
% set to 1 if in low SNR regime
ops.quadBlocks = getOr(ops, {'quadBlocks'}, 0);

orig_class = class(data);

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
    

% within frame interpolation from numBlocks -> Ly
if ops.quadBlocks    
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
    dx = xyMask(:,1:numBlocks) * squeeze(ds(:,2,:))';
    dy = xyMask(:,1:numBlocks) * squeeze(ds(:,1,:))';   
end

ds0 = cat(3, gather_try(dy), gather_try(dx));
ds0 = permute(ds0, [3 1 2]);

dx = round(dx);
dy = round(dy);

idy = repmat([1:Ly]', 1, Lx);
idx = repmat([1:Lx],  Ly, 1) ;

% shift data by smoothed shifts
dreg = zeros(size(data), orig_class);
Valid = true(Ly, Lx);
for i = 1:NT
    Im = data(:,:,i);    
    
    DX = repmat(dx(:,i),1,Lx) + idx;
    DY = repmat(dy(:,i),1,Lx) + idy;
    
    xyInvalid = DX<1 | DX>Lx | DY<1 | DY>Ly;
    Valid(xyInvalid) = false;
    
    DX(xyInvalid) = 1;
    DY(xyInvalid) = 1;
    
%     DX = mod(DX, Lx);
%     DY = mod(DY-1, Ly) + 1;
%     
    ind = DY + (DX-1) * Ly;
    Im = Im(ind);
    dreg(:,:,i) = Im;
    
end

