% takes NimgFirstRegistration mean image and aligns it to itself
% returns mean image for registration

function ops = blockAlignIterative(data, ops)

% uu = squeeze(sum(sum(data(:,:,:).^2,1),2));
% [~, isort] = sort(uu, 'descend');
% ops.mimg        = data(:,:,isort(50));


% take most correlated frames to compute initial mean image
mimg = pick_reg_init(data);

numBlocks = numel(ops.yBL);
dsold = zeros(size(data,3), 2, numBlocks);
err = zeros(ops.NiterPrealign, 1);
%%
tempSubPixel = ops.SubPixel;
for i = 1:ops.NiterPrealign    
    dsnew = zeros(size(data,3), 2, numBlocks,'double');
    Corr = zeros(size(data,3), numBlocks,'double');
    for ib = 1:numBlocks
        % collect ds
        ops.mimg = mimg(ops.yBL{ib},ops.xBL{ib});
        % compute offsets from each frame to mean image and subtract mean
        % offset (centered at offset = 0)
        [dsnew(:,:,ib), Corr(:,ib)]  = ...
            regoffKriging(data(ops.yBL{ib},ops.xBL{ib},:), ops, 1);
    end
    
    % register frames
    dreg = blockRegisterMovie(data, ops.xyMask, dsnew);
    
    % sort by correlation and take mean of most correlated frames
    [~, igood] = sort(mean(Corr,2), 'descend');
    if i<floor(ops.NiterPrealign/2)        
        igood = igood(1:100);  
    else
        igood = igood(1:round(size(data,3)/2));  
    end
    mimg = mean(dreg(:,:,igood),3);
    
    err(i) = mean(mean(sum((dsold - dsnew).^2,2),3)).^.5;
    
    dsold = dsnew;
end
ops.mimg = mimg;
ops.SubPixel = tempSubPixel;

ops.AlignNanThresh = median(Corr) - 4*std(Corr);
ops.ErrorInitialAlign = err;
ops.dsprealign = dsnew;

end 
