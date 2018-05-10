% correlate frames with each other
% sort frames by correlations 
% take mean image as mean over most correlated frames

function mimg = pick_reg_init(data)
%%
dd = bsxfun(@minus, data, mean(mean(data,1),2));

% WHITENING???
% for i = 1:size(data,3)
%    d0 = dd(:,:,i); 
%    fd0 = fft2(d0);
%    dd(:,:,i) = real(ifft2(fd0./abs(fd0)));
% end

dd = reshape(dd, [], size(dd,ndims(dd)));

CC = dd'*dd;
CC = CC./(diag(CC) * diag(CC)').^.5;

% CC = corrcoef(dd);

[CCsort, isort] = sort(CC, 2, 'descend');
nc = size(CCsort,2);
bestCC = mean(CCsort(:, 1:min(nc,20)), 2);
[~, imax] = max(bestCC);

mimg = mean(data(:,:,isort(imax, 1:min(nc,20))), 3);
end
