function [Uall,dat] = LoadMultiDay(db,iplane,iExpts, root);

ik = 0;
clear U;
nComps = 1000;
npix = 512;
clear Uall;

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
    U0 = zeros(npix,npix,nComps,'single');
    U0(dat{ik}.ops.yrange,dat{ik}.ops.xrange,:) = dat{ik}.U;
    dat{ik}.U=[];
    Uall(:,:,:,ik) = U0;
    clear U0;
end
