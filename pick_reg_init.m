
function mimg = pick_reg_init(data)

dd = data;

for i = 1:size(data,3)
   d0 = dd(:,:,i); 
   fd0 = fft2(d0);
   dd(:,:,i) = real(ifft2(fd0./abs(fd0)));
end

dd = reshape(dd, [], size(dd,3));
CC = corrcoef(dd);

[CCsort, isort] = sort(CC, 2, 'descend');

bestCC = mean(CCsort(:, 1:20), 2);
[~, imax] = max(bestCC);

mimg = mean(data(:,:,isort(imax, 1:20)), 3);
end
