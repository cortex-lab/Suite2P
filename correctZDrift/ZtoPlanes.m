% stretches Z-stack in YZ axis by angle ang
% finds best alignment of multi-plane imaging with stretched z-stack
% ang is in radians
% B is a cell array with each cell containing the mean image of that plane
function   [cMAX, ix, cZ0] = ZtoPlanes(A0, B, ang, useGPU)

[nX, nY, nZ] = size(A0);

if abs(ang) > 1e-4
    Ap = stretchZstack(A0, ang); % around axis of tissue
else
    Ap = A0;
end
% after stretching, some planes will have zeros
% exclude planes with more than 50 lines of zeros
nz = sum(squeeze(sum(Ap==0,2)==size(Ap,2)),1);
inZ = nz < 50;
Ap = Ap(:,:,inZ);

if useGPU == 1
    inZ = gather(inZ);
end
% fft of z-stack
m1 = fft(fft(Ap,[],1),[],2);
eps0 = single(1e-20);
m1 = m1./(abs(m1)+eps0);

% pad B to make it the same size as Ap
B0 = [];
for j = 1:numel(B)
    if ~isempty(B{j})
        np1 = nX-size(B{j},1); 
        np2 = nY-size(B{j},2);
        Bj = pad3Dzeros(B{j},np1,np2,0);
    
        if useGPU
            Bj = gpuArray(single(Bj));
        end
        m2 = fft(fft(Bj,[],1),[],2);
        m2 = m2./(abs(m2)+eps0);
    
        
        [c1,ix1,cz1] = ZRegPlane(m1,m2,[1:size(m1,3)],useGPU);
    
        cMAX(j) = c1;
        ix(j,:) = ix1;
        cZ(:,j) = cz1;
       
    end
end

ix(:,3)      = find(inZ,1)-1 + ix(:,3);

cZ0        = NaN*ones(size(A0,3), numel(B));
cZ0(inZ,:)   = cZ;
