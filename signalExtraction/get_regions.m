function stat = get_regions(stat, res)

Ly = res.Ly;
Lx = res.Lx;

% for computing stats
xs = repmat(1:Lx, Ly, 1);
ys = repmat((1:Ly)', 1, Lx);

xlx         = repmat(-ceil(Lx/2):1:ceil(Lx/2), 2*ceil(Lx/2)+1, 1);
rgrid       = sqrt(xlx.^2 + xlx'.^2);
rgridsort   = sort(rgrid(:), 'ascend');

% for computing neighboring pixels
T = zeros(Ly+2, Lx+2);
T(2:end-1, 2:end-1) = reshape(1:(Lx*Ly), Ly, Lx);
neigh = cat(3, T(1:end-2,2:end-1), T(3:end,2:end-1), T(2:end-1,1:end-2), T(2:end-1,3:end));
neigh = reshape(neigh, [], 4);


Nk = numel(unique(res.iclust));
for k = 1:Nk
    % needs stat, res, neigh
    pixall = stat(k).ipix;

%     minV = clustrules.parent.minPixRelVar * mean(res.M(pixall));
%     pixall(res.M(pixall)< minV) = [];
    
    whclust = 0;
    region = [];
    while numel(pixall)>0
        pixi = pixall(end);
        pixall(end) = [];
        clustn = neigh(pixi, :)';
        while 1
            newmemb = find(ismember(pixall, clustn));
            if ~isempty(newmemb)
                pixi = [pixi; pixall(newmemb)];
                newclustn = neigh(pixall(newmemb), :);
                clustn = unique(cat(1, clustn, newclustn(:)));
                pixall(newmemb) = [];
            else
                break;
            end
        end
        whclust = whclust + 1;
        
        x0 = xs(pixi); y0 = ys(pixi);
        
        rs = ((x0 - mean(x0)).^2 + (y0 - mean(y0)).^2).^.5;
        region(whclust).mrs     = mean(rs);
        region(whclust).npix    = numel(pixi);
        region(whclust).mrs0    = mean(rgridsort(1:region(whclust).npix));
        region(whclust).med     = [mean(y0) mean(x0)];
        region(whclust).ipix    = pixi;
        region(whclust).lambda  = res.lambda(pixi);
        region(whclust).V       = sum(res.M(pixi));
        region(whclust).parent      = k;  
    end
    stat(k).region = region;
    
    %     stat(k).iscell = stat(k).mrs/stat(k).mrs0<1.3 & ...
    %         stat(k).npix>50 & stat(k).npix<300;
end