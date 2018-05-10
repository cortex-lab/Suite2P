% takes NimgFirstRegistration mean image and aligns it to itself
% returns mean image for registration (ops.mimg)
function ops = alignIterative(data, ops)

fracImgPreAlign = getOr(ops, 'fracImgPreAlign', 1/2);
maxImgPreAlign = round(size(data,3) * fracImgPreAlign);

% take most correlated frames to compute initial mean image
ops.mimg = pick_reg_init(data);

dsold = zeros(size(data,3), 2);
err = zeros(ops.NiterPrealign, 1);
%%
for i = 1:ops.NiterPrealign    
  
    if ops.kriging 
        [dsnew, Corr]  = regoffKriging(data, ops, 1);
    else
        [dsnew, Corr]  = regoffLinear(data, ops, 1);
    end
    
    dreg  = rigidRegFrames(data, ops, dsnew);
    [~, igood] = sort(Corr, 'descend');
    if i<floor(ops.NiterPrealign/2)        
        igood = igood(1:min(numel(igood),100));  
    else
        igood = igood(1:maxImgPreAlign);  
    end
    ops.mimg = mean(dreg(:,:,igood),3);
    
    err(i) = mean(sum((dsold - dsnew).^2,2)).^.5;
        
    dsold = dsnew;
end

ops.AlignNanThresh = median(Corr) - 4*std(Corr);
ops.ErrorInitialAlign = err;
ops.dsprealign = dsnew;


ops.Ly = size(data,1);
ops.Lx = size(data,2);
end 
