function [Uall,dat] = LoadMultiDay(db,iplane,iExpts,root)

ik = 0;
clear U;
nComps = 1000;
clear Uall;
npix = 512;

Uall = zeros(npix,npix,nComps,length(iExpts),'single');
   
for iexp = iExpts
    CharSubDirs = '';
    for i = 1:length(db(iexp).expts)
        CharSubDirs = [CharSubDirs num2str(db(iexp).expts(i)) '_'];
    end
    blkstring = CharSubDirs(1:end-1);
    
    datexp = db(iexp).date;
    mname  = db(iexp).mouse_name;
    
    %         dat{i} = load(sprintf('~/DATA/zStacks/SVDs/%s',fs(i).name));
    fname = sprintf('SVDroi_%s_%s_plane%d.mat', mname, datexp, iplane);
    ik = ik+1;
    dat{ik} = load(fullfile(root, mname, datexp, blkstring, fname));

    % pad matrix with mirror image
    yrange = dat{ik}.ops.yrange;
    xrange = dat{ik}.ops.xrange;
    U0 = padarray(dat{ik}.U,[yrange(1)-1 xrange(1)-1],'symmetric','pre');
    U0 = padarray(U0,[npix-yrange(end) npix-xrange(end)],'symmetric','post');

    dat{ik}.ops.Nk0 = 1500;    
    dat{ik}.ops.Nk  = 1500;
    
    dat{ik}.U=[];
    Uall(:,:,:,ik) = U0;
    clear U0;
end
