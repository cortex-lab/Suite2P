function identify_redcells_sourcery(ops)

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
    [Ny Nx] = size(mimgR);
    
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
    dat.ops.mimgREDcorrected(dat.ops.yrange, dat.ops.xrange);
    %%
    
    %%%%% compute overlap with pixel map
    % (exclude pixels from other cells (cellPix))
    icell  = find(iscell);
    cellPix = false(numel(mimgR),1);
    for j = icell
        cellPix(dat.stat(j).ipix) = 1;
    end
    
    [ops, neuropMasks] = createNeuropilMasks(ops, stat, cellPix);
    
    redPix = zeros(length(dat.stat),2);
    [xx,yy] = meshgrid([1:size(mimgR0,1)],[1:size(mimgR0,2)]);
    xx=xx(:); yy=yy(:);
    xy2 = xx.^2+yy.^2;
    
    for j = 1:numel(dat.stat)
        ipix                 = dat.stat(j).ipix;
        imgpix               = false(size(mimgR0));
        imgpix(ipix)         = 1;
        rpix                 = mimgR0(ipix);%.*(dat.stat(j).lambda'/sum(dat.stat(j).lambda));
        [ix,iy]              = ind2sub(size(mimgR0),ipix);
        params               = FitCircle(ix,iy);
        cellradius           = params.ra;
        yc                   = params.xc;
        xc                   = params.yc;
        icircle              = (xy2 - 2*xc*xx - 2*yc*yy + xc^2 + yc^2) < (cellradius+15)^2;
        extpix               = icircle & ~imgpix(:) & ~cellPix;
        ext_rpix             = mimgR0(extpix);%/length(extpix);
        redPix(j,:)          = [mean(rpix) mean(ext_rpix)];
    end
    redPix = redPix - min(redPix(:));
    %     redpix(~dat.cl.isroi,:) = NaN;
    %%

    % set threshold for redpix
    if isfield(ops,'redthres')
        redthres = ops.redthres;
    else
        redthres = 1.5;
    end
    if isfield(ops,'redmax')
        redmax = ops.redmax;
    else
        redmax = 1;
    end
    rrat = redPix(:,1)./(redPix(:,2)+redPix(:,1));
    redcell  = rrat > nanmean(rrat)+redthres*nanstd(rrat);
    notred   = rrat <= nanmean(rrat) + redmax*nanstd(rrat);
    
    fprintf('plane %d  reds %d\n',iplane,sum(redcell(:)&iscell(:)));
    
%     dat.cl.redcell = redcell(:);
%     dat.cl.notred  = notred(:);

    for j = 1:length(dat.stat)
        dat.stat(j).redcell = redcell(j);
        dat.stat(j).redprob = rrat(j);
    end

    save(fname, '-struct', 'dat')
end









