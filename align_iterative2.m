function ops = align_iterative2(data, ops)

% uu = squeeze(sum(sum(data(:,:,:).^2,1),2));
% [~, isort] = sort(uu, 'descend');
% ops.mimg        = data(:,:,isort(50));


ops.mimg = pick_reg_init(data);

dsold = zeros(size(data,3), 2);
err = zeros(ops.NiterPrealign, 1);
%%
ops.SubPixel = Inf;
for i = 1:ops.NiterPrealign    
    
    [dsnew, Corr]  = registration_offsets(data, ops, 1);
    dreg  = register_movie(data, ops, dsnew);
    [~, igood] = sort(Corr, 'descend');
    
    ops.mimg = mean(dreg(:,:,igood),3);
    
    err(i) = mean(sum((dsold - dsnew).^2,2)).^.5;
        
    dsold = dsnew;
end

ops.AlignNanThresh = median(Corr) - 4*std(Corr);
ops.ErrorInitialAlign = err;
ops.dsprealign = dsnew;

end 
