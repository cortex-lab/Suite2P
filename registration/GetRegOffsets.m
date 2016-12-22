function [dsall, ops1] = GetRegOffsets(data, j, iplane0, ops, ops1, yFOVs, xFOVs)

dsall = zeros(size(data,3), 2, size(xFOVs,2));
for i = 1:numel(ops.planesToProcess)
    ifr0 = iplane0(ops.planesToProcess(i));
    indframes = ifr0:ops.nplanes:size(data,3);
    
    for l = 1:size(xFOVs,2)
        dat = data(yFOVs(:,l),xFOVs(:,l),indframes);
        if ~isempty(ops.smooth_time_space)
            dat = smooth_movie(dat, ops);
        end
        [ds, Corr]  = regoffKriging(dat, ops1{i,l}, 0);
        
        %ds          = RemoveBadShifts(ds);
        
        dsall(indframes,:, l)  = ds;
        % collect ds
        if j==1
            ds(1,:,:) = 0;
        end
        ops1{i,l}.DS          = cat(1, ops1{i,l}.DS, ds);
        ops1{i,l}.CorrFrame   = cat(1, ops1{i,l}.CorrFrame, Corr);
    end
    
    1;
    % check if there was a sharp drop in fluorescence
%     lbright = sq(mean(data(:,:,indframes),2));
%     mlbright = mean(ops1{i}.mimg, 2);
%     
%     lbright = bsxfun(@rdivide, lbright, mlbright);
%     badi = max(abs(lbright(1:end-4,:) - lbright(5:end,:)), [], 1) > .5;
%     badi = find(badi);
%     
%     ops1{i}.badframes(sum(ops1{i}.Nframes) + badi) = true;    
end