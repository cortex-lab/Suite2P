function dat = LoadRegops(db,iplane,iExpts,root)

ik = 0;
clear U;
nComps = 1000;
clear Uall;
npix = 512;
   
for iexp = iExpts
    CharSubDirs = '';
    for i = 1:length(db(iexp).expts)
        CharSubDirs = [CharSubDirs num2str(db(iexp).expts(i)) '_'];
    end
    blkstring = CharSubDirs(1:end-1);
    
    datexp = db(iexp).date;
    mname  = db(iexp).mouse_name;
    
    %         dat{i} = load(sprintf('~/DATA/zStacks/SVDs/%s',fs(i).name));
    fname = sprintf('regops_%s_%s_plane%d.mat', mname, datexp, iplane);
    ik = ik+1;
    dat{ik} = load(fullfile(root, mname, datexp, blkstring, fname));

end
