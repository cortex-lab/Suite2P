% function run_pipeline(db, ops0, clustrules)

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
for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    ops     = ops1{i};

    ops.iplane  = iplane;

    if numel(ops.yrange)>ops.Ly/2 && numel(ops.xrange)>ops.Lx/2
        if getOr(ops, {'getSVDcomps'}, 0)
            ops    = get_svdcomps(ops);
        end
        if ops.getROIs || getOr(ops, {'writeSVDroi'}, 0)
            [ops, U, Sv]        = get_svdForROI(ops);
        end
        if ops.getROIs
            U =  reshape(U, [], size(U,3));
            
            switch clustModel
                case 'standard'
                    [ops, stat, res]  = fast_clustering(ops,U, Sv);
                case 'neuropil'                    
                    [ops, stat, res]  = fast_clustering_with_neuropil(ops,U, Sv);
            end
            
            apply_ROIrules(ops, stat, res, clustrules);

            switch neuropilSub
                case 'surround'
                    get_signals_and_neuropil(ops, iplane);
                case 'none'
                    get_signals(ops, iplane);
                case 'model'
                     get_signals_NEUmodel(ops, iplane);
            end
        end
    end

    if ops.DeleteBin
        fclose('all');
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;