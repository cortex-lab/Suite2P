clear all;
% path to db files
addpath('/media/carsen/DATA1/2P/dbfiles/')
compile_dbs;

% path to saved F files (and root directory)
ops0.RootDir                = '//zserver.cortexlab.net/Data/Subjects/';
ops0.ResultsSavePath        = '/media/carsen/DATA1/2P/F/';

% choose experiment
iexp = 31;
db = db0(iexp);

% build ops
ops = build_ops3(db, ops0);
root = ops.ResultsSavePath;


% load mean image and stats for each plane in recording
stat = [];
xL = [];
yL = [];
clear mimg;
for iplane = 1:ops.nplanes
    fname  = sprintf('F_%s_%s_plane%d.mat', db.mouse_name, db.date, iplane);
    dat = load(fullfile(root, fname));

    mimg{iplane} = dat.ops.mimg1(dat.ops.yrange, dat.ops.xrange);
        
    stat = cat(2, stat, dat.stat);
    yL(iplane) = numel(dat.ops.yrange);
    xL(iplane) = numel(dat.ops.xrange);
end


% compute outlines of cells
clear img;
ij = 1;
clf;
for iplane = 1:ops.nplanes
    mimg{iplane} = mimg{iplane} - min(mimg{iplane}(:));
    mimg{iplane} = mimg{iplane} / max(mimg{iplane}(:));
    
    % background of mean image
    img{iplane} = repmat(mimg{iplane}*2,1,1,3);
    img{iplane} = rgb2hsv(img{iplane});
    img{iplane} = reshape(img{iplane},[],3);
    
    colormap('gray');
    hold all;
    while stat(ij).iplane == iplane 
        if stat(ij).iscell
            % find pixels that are exterior to the cell 
            idist  = sqrt(bsxfun(@minus, stat(ij).xpix', stat(ij).xpix).^2 + ...
                bsxfun(@minus, stat(ij).ypix', stat(ij).ypix).^2);
            idist  = idist - diag(NaN*diag(idist));
            extpix = sum(idist <= sqrt(2)) <= 6;
            xext = stat(ij).xpix(extpix);
            yext = stat(ij).ypix(extpix);
                
            % color the exterior pixels
            ipix = sub2ind([yL(iplane) xL(iplane)], yext, xext);
            img{iplane}(ipix, 1) = rand;
            img{iplane}(ipix, 2) = 1;
            img{iplane}(ipix, 3) = 1;
            
            
        end
        
        ij = ij+1;
        
        if ij > numel(stat)
            break;
        end
    end
    
end

%%
clf;
set(gcf,'color','w');

% plot of one plane with circled ROIs
iplane = 5;

imgj = img{iplane};
imagesc(hsv2rgb(reshape(imgj,yL(iplane),xL(iplane),3)))
axis off;
axis square;



    
    
    
    
  
