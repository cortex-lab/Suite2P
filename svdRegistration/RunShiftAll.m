function [Uall,dat,ops,refShift] = RunShiftAll(Uall,dat,nBest,isGPU)



isSVD = 1;
npix = size(Uall,1);
pixShifts = zeros(npix,npix,2,length(dat));
nDat = length(dat);
nComps = 10;        %%%% <--- how many svd components to use for estimates


% Uall = Uall(:,:,1:nComps,:);
% Ufft = fft(fft(Uall, [], 1), [], 2);
% Ufft = Ufft./(abs(Ufft) + 1e-4);
% Uall = real(ifft(ifft(Ufft, [], 1), [], 2));

%% initialize shifts with best dataset for all
yrange0 = dat{nBest}.ops.yrange;
xrange0 = dat{nBest}.ops.xrange;
dat{nBest}.ops.yrange2=yrange0;
dat{nBest}.ops.xrange2=xrange0;
dinds = true(1,nDat);
dinds(nBest) = 0;
dinds = find(dinds);
% reference image
A = single(Uall(:,:,1:nComps,nBest));

nIter = 1;
for i = 1:nIter
    for nD = 1:nDat
        B = single(Uall(:,:,1:nComps,nD));
        % compute quadrant based shifts
        [ops,pixShift,~,~,cX] = QuadrantPixelShift(dat{nD}.ops, ...
            npix,A,B,isGPU,isSVD);
        pixShifts(:,:,:,nD) = pixShift;
    end
    % normalize pixShifts so middle of image is constant
    pixShiftsY = squeeze(pixShifts(:,:,1,:));
    pixShiftsX = squeeze(pixShifts(:,:,2,:));
    pixShiftsY = mean(pixShiftsY(:));
    pixShiftsX = mean(pixShiftsX(:));
    pixShifts(:,:,1,:) = pixShifts(:,:,1,:) - pixShiftsY;
    pixShifts(:,:,2,:) = pixShifts(:,:,2,:) - pixShiftsX;
    
    pixShiftMean = ones(npix,npix,2,'single');
    pixShiftMean(:,:,1) = -1*pixShiftsY;
    pixShiftMean(:,:,2) = -1*pixShiftsX;
    
    if nIter>1
        % shift SVDs
        CX = [];
        A = ShiftSVD(A, pixShiftMean,npix,isGPU);
        Un = zeros(npix,npix,nComps,'single');
        for nD = 1:nDat
            U0 = ShiftSVD(Uall(:,:,1:nComps,nD), pixShifts(:,:,:,nD),npix,isGPU);
            Uall(:,:,1:nComps,nD) = U0;
            [U0,cX] = SVDProj(A,U0(:,:,1:nComps),isGPU);
            Un = Un + U0;
            CX(nD) = cX;
        end
        disp([mean(CX) CX])
        Un = Un/nDat;
        Un2 = reshape(Un,[],nComps);
        [u,s,v] = svd(Un2'*Un2);
        Un2 = normc(Un2*u);
        Un  = reshape(Un2,npix,npix,nComps);
        
        A = Un;
        
        imagesc(pixShifts(:,:,1,1))
        colorbar
        drawnow
    end
    %clear Un;
end
%%

for nD=1:nDat
    fprintf('dat %d, yshift %2.2f xshift %2.2f\n',...
        nD,mean(mean(pixShifts(:,:,1,nD))),mean(mean(pixShifts(:,:,2,nD))))
end
refShift = A;
clear Un;



% shift U's
for nD = 1:nDat
%     dat{nD}.ops = ops;
    pixShift = pixShifts(:,:,:,nD);
    dat{nD}.pixShift = pixShift;
    
    % shift B to A
    % aligned new U to reference image
    yrange = dat{nD}.ops.yrange;
    xrange = dat{nD}.ops.xrange;
    ymax = round(max(0,max(max(pixShift(:,:,1)))));
    xmax = round(max(0,max(max(pixShift(:,:,2)))));
    ymin = round(min(0,min(min(pixShift(:,:,1)))));
    xmin = round(min(0,min(min(pixShift(:,:,2)))));
    dat{nD}.ops.yrange2 = [yrange(1)+ymax : yrange(end)+ymin];
    dat{nD}.ops.xrange2 = [xrange(1)+xmax : xrange(end)+xmin];
    
    U0 = ShiftSVD(Uall(:,:,:,nD), pixShift,npix,isGPU);
    Uall(:,:,:,nD) = U0;
    clear U0;
end

% intersect all indices
yrange = dat{1}.ops.yrange2;
xrange = dat{1}.ops.xrange2;

for nD = 2:length(dat)
    yrange = intersect(yrange,dat{nD}.ops.yrange2);
    xrange = intersect(xrange,dat{nD}.ops.xrange2);
end
ops.yrange = yrange;
ops.xrange = xrange;

for nD = [1:nDat]
    dat{nD}.ops.yrange0 =  yrange;
    dat{nD}.ops.xrange0 =  xrange;
end

for nD = 1:nDat
    for i = 1:size(Uall,3)
        Uall(:,:,i,nD) = Uall(:,:,i,nD) * dat{nD}.Sv(i)^.5;
    end
end
