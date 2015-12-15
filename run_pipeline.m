function run_pipeline(iexp, db, ops0, clustrules)

if ops0.CopyDataLocally
    db0 = copy_from_zserver(db(iexp), ops0);
    ops = build_ops2(db0, ops0);
else
    ops = build_ops(db(iexp), ops0);
end

if ops.useGPU
    gpuDevice(1);   % reset GPU at each dataset
end
%
ops1         = reg2Pnew(ops);  % do registration

for i = 1:length(ops.planesToProcess)
    iplane  = ops.planesToProcess(i);
    ops     = ops1{i};
    ops.iplane  = iplane;
    
    if numel(ops.yrange)>300 && numel(ops.xrange)>300
        if ops.getSVDcomps
            ops    = get_svdcomps(ops);
        end
        
        if ops.getROIs            
            [ops, U, Sv]        = get_svdForROI(ops);
            %
            [ops, stat0, res0]  = fast_clustering(ops, reshape(U, [], size(U,3)), Sv);
            %
            apply_ROIrules(ops, stat0, res0, clustrules);
            %
            get_signals_and_neuropil(ops, iplane);
        end
    end
    
    if ops.DeleteBin
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;