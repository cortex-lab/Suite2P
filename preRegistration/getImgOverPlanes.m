function newTargets = getImgOverPlanes(ops1)
nFOV = size(ops1,2);
Ly = ops1{1}.Ly;
Lx = ops1{1}.Lx;
planesToInterpolate = ops1{1}.planesToInterpolate;
ds = zeros(length(planesToInterpolate), 2, nFOV);
% planesToInterpolate should be indices of non-flyback planes
% (only those should be aligned)
for i = 2:length(planesToInterpolate) % align consecutive planes
    pl = planesToInterpolate(i);
    for j = 1:nFOV
        ds(i,:,j) = registration_offsets(ops1{pl,j}.mimg, ops1{pl-1,j}, 0);
    end
end
ds = cumsum(ds,1);
ds = bsxfun(@minus, ds, mean(ds,1));
images = zeros(Ly, Lx, length(planesToInterpolate), nFOV);
for i = 1:length(planesToInterpolate)
    for j = 1:nFOV
        images(:,:,i,j) = ops1{planesToInterpolate(i),j}.mimg;
    end
end
ds = reshape(permute(ds, [1 3 2]), [], 2);
images = reshape(images, Ly, Lx, []);
newTargets = register_movie(images, ops1{1}, ds);
newTargets = reshape(newTargets, Ly, Lx, length(planesToInterpolate), nFOV);