function identify_redcells_sourcery(db, ops0)

ops = build_ops3(db, ops0);

redcells = [];

for i = 1:length(ops.planesToProcess)
    
    iplane  = ops.planesToProcess(i);
    %%
    try
        fname = sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane);
        dat = load(fname);
        while isfield(dat, 'dat')
            dat = dat.dat;
        end
    catch
        fname = sprintf('%s/F_%s_%s_plane%d_proc.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        while isfield(dat, 'dat')
            dat = dat.dat;
        end
    end
    
    % subtract bleedthrough of green into red channel
    % if there's a specific green channel (from a red+green short session)
    % using non-rigid linear regression
    if isfield(dat.ops,'mimgGREEN')
        mimgG = dat.ops.mimgGREEN;
    else
        mimgG = dat.ops.mimg1;
    end
    mimgR = dat.ops.mimgRED;
    mimgR = mimgR(dat.ops.yrange, dat.ops.xrange);
    mimgG = mimgG(dat.ops.yrange, dat.ops.xrange);
    [Ny, Nx] = size(mimgR);
    
    % regression
    nblks               = 3;
    yB                  = round(linspace(1,Ny,nblks+1));
    xB                  = round(linspace(1,Nx,nblks+1));
    [yBL,xBL] = MakeQuadrants(yB,xB);
    msk                 = zeros(Ny,Nx,length(yB));
    for j = 1:length(yBL)
        xg              = mimgG(yBL{j},xBL{j});
        xr              = mimgR(yBL{j},xBL{j});
        A               = polyfit(xg,xr,1);
        msk(:,:,j)      = QuadrantMask(yBL,xBL,Ny,Nx,j);
        gw(j)           = A(1);
    end
    msk                 = bsxfun(@times,msk,1./sum(msk,3));
    pixweight           = zeros(Ny,Nx);
    for j = 1:length(yBL)
        pixweight = pixweight + msk(:,:,j)*gw(j);
    end
    mimgR0              = mimgR - pixweight.*mimgG;
    
    % save to dat for GUI
    dat.ops.mimgREDcorrected = zeros(dat.ops.Ly, dat.ops.Lx);
    dat.ops.mimgREDcorrected(dat.ops.yrange, dat.ops.xrange) = mimgR0;
    %
    
    %%
    %%%% compute red pixels in cell and in area around cell
    % (exclude pixels from other cells (cellPix))
    % (use surround neuropil masks to compute red around cell)
    ops.ratioNeuropil     = 4;
    ops.minNeuropilPixels = 80;
    [dat.stat, cellPix, ~]       = createCellMasks(dat.stat, Ny, Nx);
    [~, neuropMasks] = createNeuropilMasks(ops, dat.stat, cellPix);
    
    %
    redSum = zeros(length(dat.stat),2);
    for k = 1:numel(dat.stat)
        ipix                 = dat.stat(k).ipix;
        rpix                 = mimgR0(ipix) .* dat.stat(k).lam/sum(dat.stat(k).lam);
        % rescale by lam? .*(dat.stat(j).lambda'/sum(dat.stat(j).lambda));
        extpix               = squeeze(neuropMasks(k,:,:));
        extpix               = find(extpix~=0);
        ext_rpix             = mimgR0(extpix);
        redSum(k,:)          = [sum(rpix) mean(ext_rpix)];
    end
    redSum = redSum - min(redSum(:));
    
    % set threshold for redpix
    ops.redthres = getOr(ops, 'redthres', 1.5);
    ops.redmax   = getOr(ops, 'redmax', 1.25);
    
    rrat = redSum(:,1)./(redSum(:,2)+redSum(:,1));
    redcell  = rrat > nanmean(rrat) + ops.redthres*nanstd(rrat);
    notred   = rrat <= nanmean(rrat) + ops.redmax*nanstd(rrat);
  
    
    fprintf('plane %d  reds %d\n',iplane,sum(redcell(:)));
    

    %%
    for j = 1:length(dat.stat)
        dat.stat(j).redcell = redcell(j);
        dat.stat(j).redprob = rrat(j);
        dat.stat(j).notred  = notred(j);
    end

    dat.ops.redthres = ops.redthres;
    dat.ops.redmax   = ops.redmax;
    save(fname, '-struct', 'dat')
end


  
    
%%
if 0
    [~,ix] = sort(rrat, 'descend');
    for k = ix(:)'
        clf;
        cent = round(dat.stat(k).med);
        yl = cent(1) + [-20:20];
        yl(yl<1 | yl > Ny) = [];
        xl = cent(2) + [-20:20];
        xl(xl<1 | xl > Nx) = [];
        imagesc(mimgR0)
        hold all;
        plot(cent(2),cent(1),'k*');
        axis([xl(1) xl(end) yl(1) yl(end)]);
        title([rrat(k) redcell(k)]);
        drawnow;
        pause;
    end
end
    %






