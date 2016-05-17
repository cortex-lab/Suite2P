function [U,cX] = SVDProj(A,B,isGPU)

if isGPU == 1
    A = gpuArray(A);
    B = gpuArray(B);
    U = gpuArray.zeros(size(A),'single');
else
    U = zeros(size(A),'single');
end

[nX nY nComps] = size(A);
for i = 1:nComps
    for j = 1:nComps
        U(:,:,i) = U(:,:,i) + B(:,:,j) * (sum(sum(B(:,:,j) .* A(:,:,i))));
    end
end


% compare U and A
cX = 0;
for i = 1:nComps
    for j = 1:nComps
%         [cX0,ix,cc] = regZ(U(:,:,i),A(:,:,j));
         cX0 = sum(sum((U(:,:,i).*A(:,:,j))));
         cX = cX + cX0;
    end
end
cX = cX/nComps;


if isGPU == 1
	cX = gather(cX);
	U = gather(U);
end
