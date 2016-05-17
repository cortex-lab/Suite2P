function [dat] = ShiftAllMasks(res,ops,isGPU);

B = single(dat{nBest}.ops.mimg1);
for nD = dinds
    A = single(dat{nD}.ops.mimg1);
    [res0,pixInv] = ShiftMasks(res,ops,A,B,isGPU);
    dat{nD}.res    = res0;
    dat{nD}.pixInv = pixInv;
end
