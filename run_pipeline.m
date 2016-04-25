function run_pipeline(db, ops0, clustrules)

ops = build_ops3(db, ops0);

if ops.useGPU
    gpuDevice(1);   % reset GPU at each dataset
end
%

clustModel     = getOr(ops, {'clustModel'}, 'standard');
neuropilSub    = getOr(ops, {'neuropilSub'}, 'surround');
splitBlocks    = getOr(ops, {'splitBlocks'}, 'none');

if iscell(splitBlocks)
    ops1         = blockReg2P(ops);  % do registration
else
    ops1         = reg2P(ops);  % do registration
end
 %%
for i = 3:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    ops     = ops1{i};
    
    ops.iplane  = iplane;
    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
        continue;
    end
    ops.Nk = 1500;
    ops.Nk0 = 1500;
    
    if getOr(ops, {'getSVDcomps'}, 0)
        ops    = get_svdcomps(ops);
    end
    if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
        [ops, U, Sv]    = get_svdForROI(ops);
        ops.U           = ops.U(:,:,1:20);
        ops.Sv          = ops.Sv(1:20);
    end
    if ops.getROIs
        switch clustModel
            case 'standard'
                [ops, stat, res]  = fast_clustering(ops,U, Sv);
            case 'neuropil'
                [ops, stat, res]  = fast_clustering_with_neuropil(ops,U, Sv);
        end
        
        [stat, res] = apply_ROIrules(ops, stat, res, clustrules);
        
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