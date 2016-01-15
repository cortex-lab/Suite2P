function run_pipeline(db, ops0, clustrules)

ops = build_ops3(db, ops0);

if ops.useGPU
    gpuDevice(1);   % reset GPU at each dataset
end
%
ops1         = reg2P(ops);  % do registration

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
<<<<<<< HEAD
%             get_signals_and_neuropil(ops, iplane);
            get_signals(ops, iplane);
=======
            %get_signals_and_neuropil(ops, iplane);
			get_signals_and_neuropil(ops, iplane);
>>>>>>> origin/master
        end
    end
    
    if ops.DeleteBin
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose all;