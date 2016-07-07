function  run_pipeline(db, ops0, clustrules)


% ops0.TileFactor (or db(iexp).TileFactor) can be set to multiply the number of default tiles for the neuropil

ops0.nimgbegend                     = getOr(ops0, {'nimgbegend'}, 250);
ops0.LoadRegMean                    = getOr(ops0, {'LoadRegMean'}, 0);
ops0.NiterPrealign                  = getOr(ops0, {'NiterPrealign'}, 10);
ops0.registrationUpsample           = getOr(ops0, {'registrationUpsample'}, 1);  % upsampling factor during registration, 1 for no upsampling is much faster, 2 may give better subpixel accuracy
ops0.niterclustering                = getOr(ops0, {'niterclustering'}, 50);   % how many iterations of clustering
ops0.getROIs                        = getOr(ops0, {'getROIs'}, 1);   % whether to run the optimization
ops0.getSVDcomps                    = getOr(ops0, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing
ops0.nSVD                           = getOr(ops0, {'nSVD'}, 1000);   % how many SVD components to save to disk

clustrules.npix_fraclow             = getOr(clustrules, {'npix_fraclow'}, 1/4);
clustrules.npix_frachigh            = getOr(clustrules, {'npix_frachigh'}, 20);
clustrules.diameter                 = getOr(clustrules, {'diameter'}, 10);
ops0.diameter                       = clustrules.diameter;
clustrules.MinNpix                  = round(pi/4 * clustrules.diameter^2 *clustrules.npix_fraclow);
clustrules.MaxNpix                  = round(pi/4 * clustrules.diameter^2 * clustrules.npix_frachigh);
clustrules.Compact                  = getOr(clustrules, {'Compact'}, 2);
clustrules.parent                   = getOr(clustrules, {'parent'}, []);
clustrules.parent.minPixRelVar      = getOr(clustrules.parent, {'minPixRelVar'}, 1/10);
clustrules.parent.PixelFractionThreshold     = getOr(clustrules.parent, {'PixelFractionThreshold'}, 0.5);
clustrules.parent.MaxRegions        = getOr(clustrules.parent, {'MaxRegions'}, 10);
 
ops = build_ops3(db, ops0);

if ops.useGPU
    gpuDevice(1);   % reset GPU at each dataset
end
%

clustModel     = getOr(ops, {'clustModel'}, 'standard');
neuropilSub    = getOr(ops, {'neuropilSub'}, 'surround');
splitBlocks    = getOr(ops, {'splitBlocks'}, 'none');
%
if iscell(splitBlocks)
    ops1         = blockReg2P(ops);  % do registration
else
    ops1         = reg2P(ops);  % do registration
end
%
%%
for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    ops     = ops1{i};
    %
    
    ops.iplane  = iplane;
    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
%         continue;
    end
    
    if getOr(ops, {'getSVDcomps'}, 0)
        ops    = get_svdcomps(ops);
    end
    if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
        [ops, U, Sv]    = get_svdForROI(ops);
    end
    
    if ops.getROIs
        switch clustModel
            case 'standard'
                [ops, stat, res]  = fast_clustering(ops,U, Sv);
            case 'neuropil'                    
%                 [ops, stat, res]  = fast_clustering_with_neuropil(ops,U, Sv);
                  % better model of the neuropil
                  [ops, stat, res]  = fastClustNeuropilCoef(ops,U, Sv);
        end
        
        [stat2, res2] = apply_ROIrules(ops, stat, res, clustrules);
        
        switch neuropilSub
            case 'surround'
                get_signals_and_neuropil(ops, iplane);
            case 'none'
                get_signals(ops, iplane);
            case 'model'
                get_signals_NEUmodel(ops, iplane);
        end
    end
    

    if ops.DeleteBin
        fclose('all');
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;