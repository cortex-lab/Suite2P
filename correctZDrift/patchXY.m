% make yx patches with centers at (yc, xc) and size (pixy, pixx)
% returns pixel indices x number of patches matrix    
function ipix0 = patchXY(Ly, Lx, yc, xc, pixy, pixx)
    
ipixy  = bsxfun(@plus, yc(:) - floor(pixy/2), [0:pixy-1]);
ipixx  = bsxfun(@plus, xc(:) - floor(pixx/2), [0:pixx-1]);
ipix0  = zeros(pixy*pixx, numel(yc));
for j = 1:numel(yc)
    [iy, ix] = ndgrid(ipixy(j,:), ipixx(j,:));
    ipix0(:,j)   = sub2ind([Ly Lx], iy(:), ix(:));
end
    