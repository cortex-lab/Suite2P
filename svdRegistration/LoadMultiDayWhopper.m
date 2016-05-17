function [Uall,dat] = LoadMultiDayWhopper

ik = 0;
clear U;
nComps = 1000;
npix = 512;

clear Uall;
clear dat;
fs = dir('~/DATA/zStacks/SVDs/bad_align/*.mat');
Uall = zeros(npix,npix,nComps,length(fs),'single');

for i = [1:length(fs)]
    dat{i} = load(sprintf('~/DATA/zStacks/SVDs/bad_align/%s', ...
                          fs(i).name));
    %U0 = zeros(npix,npix,nComps,'single');
    %U0(dat{i}.ops.yrange,dat{i}.ops.xrange,:) = dat{i}.U;
    yrange = dat{i}.ops.yrange;
    xrange = dat{i}.ops.xrange;

    U0 = padarray(dat{i}.U,[yrange(1)-1 xrange(1)-1],'symmetric','pre');
    U0 = padarray(U0,[npix-yrange(end) npix-xrange(end)],'symmetric','post');

    dat{i}.U=[];
    Uall(:,:,:,i) = U0;
    clear U0;
    
end

    