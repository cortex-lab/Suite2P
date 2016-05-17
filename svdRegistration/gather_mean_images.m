% master_fileMP
addpath('C:\CODE\GitHub\Suite2P\svdRegistration')

for iplane = 1:10 % [3:9]
%     clearvars -except iplane  clustrules
    isSVD = 1;
    isGPU = 1;
    nBest = 4;

    root = 'E:\DATA';
    make_db_20plane_1x;
    iExpts = [17:22];
    
    dat = LoadRegops(db,iplane,iExpts, root);
    for j = 1:length(dat)
       mig{iplane, j} = dat{j}.ops.mimg1; 
    end
end



