

for nD = 1:length(dat)
    A = single(dat{nD}.ops.mimg1);
    for nD2 = 1:length(dat)
        B = single(dat{nD2}.ops.mimg1);
        [cx,ix]=regZ(A,B);
        cxAll(nD,nD2) = cx;
    end
end