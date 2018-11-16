
function [stat, res] = trimPixels(stat, res)

Nk = numel(stat);

stat(Nk+1).mrs = Inf;
stat(Nk+1).mrs0 = 1;


for k = 1:Nk
    % trim everything below 5 pixels
   ibadreg = [stat(k).region.npix] < 10;
    
   ibad = [];
   for j = 1:length(stat(k).region)
       if ibadreg(j)
        ibad = cat(1, ibad, stat(k).region(j).ipix);
       end
   end
   stat(k).region = stat(k).region(~ibadreg);
   
   ibadpixels = ismember(stat(k).ipix, ibad);
   stat(k).ipix = stat(k).ipix(~ibadpixels);
   stat(k).lambda = stat(k).lambda(~ibadpixels);
   
   res.iclust(ibad) = Nk+1;
   
   stat(Nk+1).ipix   = cat(1, stat(Nk+1).ipix, ibad);   
end

stat(Nk+1).lambda = res.lambda(stat(Nk+1).ipix);
stat(Nk+1).npix = numel(stat(Nk+1).ipix);
    
end