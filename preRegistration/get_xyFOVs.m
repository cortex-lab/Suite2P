function [xFOVs, yFOVs] = get_xyFOVs(ops)


Ly = ops.Ly;
Lx = ops.Lx;
iR = ops.splitFOV;

ny = floor(Ly/iR(1));
nx = floor(Lx/iR(2));

xFOVs = zeros(nx, iR(2), iR(1));
yFOVs = zeros(nx, iR(2), iR(1));
for i = 1:iR(1)
    for j = 1:iR(2)
        xFOVs(:,j,i) = [1:nx] + (j-1)*nx;
        yFOVs(:,j,i) = [1:ny] + (i-1)*ny;
    end
end

xFOVs = xFOVs(:,:);
yFOVs = yFOVs(:,:);
