

ops = build_ops3(db0(iexp), ops0);
redcells = [];
for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    try
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d_proc.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        while isfield('dat', 'dat')
            dat = dat.dat;
        end
    catch
        fname = sprintf('%s/F_%s_%s_plane%d_Nk%d.mat', ops.ResultsSavePath, ops.mouse_name, ops.date, iplane, ops.Nk);
        dat = load(fname);
        while isfield('dat', 'dat')
            dat = dat.dat;
        end
    end
    
    % subtract bleedthrough of green into red channel
    % if there's a specific green channel (from a red+green short session)
    if isfield(dat.ops,'mimgGREEN')
        mimgG = dat.ops.mimgGREEN;
    else
        mimgG = dat.ops.mimg;
    end
    mimgR = dat.ops.mimgRED;
    
    % subsample mimgR and mimgG
    mimgR = mimgR(dat.ops.yrange,dat.ops.xrange);
    mimgG = mimgG(dat.ops.yrange,dat.ops.xrange);
    [nY nX] = size(mimgR);
    %
    % pixels with ROIs
    iscell = [dat.stat.mrs]./[dat.stat.mrs0]<1.2 & [dat.stat.npix]<150 & [dat.stat.npix]>20;
    icell  = find(iscell);
    cellpix = false(numel(mimgR),1);
    for j = icell
        cellpix(dat.stat(j).ipix) = 1;
    end
    
    nonrig = 1;
    if nonrig
        % non-rigid linear regression
        nblks               = 3;
        yB                  = round(linspace(1,nY,nblks+1));
        xB                  = round(linspace(1,nX,nblks+1));
        [yBL,xBL,numBlocks] = MakeQuadrants(yB,xB);
        msk                 = zeros(nY,nX,length(yB));
        for j = 1:length(yBL)
            xg              = mimgG(yBL{j},xBL{j});
            xr              = mimgR(yBL{j},xBL{j});
            A               = polyfit(xg,xr,1);
            msk(:,:,j)      = QuadrantMask(yBL,xBL,nY,nX,j);
            gw(j)           = A(1);
        end
        msk                 = bsxfun(@times,msk,1./sum(msk,3));
        pixweight           = zeros(nY,nX);
        for j = 1:length(yBL)
            pixweight = pixweight + msk(:,:,j)*gw(j);
        end
        mimgR0              = mimgR - pixweight.*mimgG;
    else
        % fit to non-ROI pixels
        mimgGn              = mimgG;
        mimgGn(cellpix)     = 0;
        mimgRn              = mimgR;
        mimgRn(cellpix)     = 0;
        
        xg                  = mimgGn;%my_conv2(mimgGn,[1 1],[1 2]);
        xr                  = mimgRn;%my_conv2(mimgRn,[1 1],[1 2]);
        A                   = polyfit(xg(:),xr(:),1);
        plot(xg(1:100:end)*A(1)+A(2),xr(1:100:end),'.')
        mimgR0              = mimgR - (A(1)*mimgG);
    end
    dat.ops.mimgREDcorrected = mimgR0;
    
    mimgR0 = normalize_image(mimgR0);
    %%%%% compute overlap with pixel map
    % (exclude pixels from other cells (cellpix))
    if ~isfield(dat.cl,'redcell')
        dat.cl.redcell = zeros(length(dat.stat), 1);
    end
    redpix = zeros(length(dat.stat),2);
    [xx,yy] = meshgrid([1:size(mimgR0,1)],[1:size(mimgR0,2)]);
    xx=xx(:); yy=yy(:);
    xy2 = xx.^2+yy.^2;
    
    for j = find(dat.cl.isroi)
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
        extpix               = icircle & ~imgpix(:) & ~cellpix;
        ext_rpix             = mimgR0(extpix);%/length(extpix);
        redpix(j,:)          = [mean(rpix) mean(ext_rpix)];
    end
    redpix = redpix - min(redpix(:));
    redpix(~dat.cl.isroi,:) = NaN;
    
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
    rrat = redpix(:,1)./(redpix(:,2)+redpix(:,1));
    redcell  = rrat > nanmean(rrat)+redthres*nanstd(rrat);
    notred   = rrat <= nanmean(rrat) + redmax*nanstd(rrat);
    
    fprintf('plane %d  reds %d\n',iplane,sum(redcell(:)&iscell(:)));
    
    dat.cl.redcell = redcell(:);
    dat.cl.notred  = notred(:);
    
    save(fname, '-struct', 'dat')
end









