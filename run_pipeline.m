function run_pipeline_new(sm)

% options.TileFactor (or db(iexp).TileFactor) can be set to multiply the number of default tiles for the neuropil

options = sm.options;
options.splitROIs                      = getOr(options, {'splitROIs'}, 1);
options.LoadRegMean                    = getOr(options, {'LoadRegMean'}, 0);
options.getROIs                        = getOr(options, {'getROIs'}, 1);   % whether to run the optimization
options.getSVDcomps                    = getOr(options, {'getSVDcomps'}, 0);   % whether to save SVD components to disk for later processing

%options                                = build_ops3(db, options);
% options initialization is done elsewhere
options.diameter = sm.getOr('diameter', 8*sm.getOr('zoom'));

if isempty(options.diameter)
    warning('you have not specified mean diameter of your ROIs')
    warning('for best performance, please set db(iexp).diameter for each experiment')
end
%options.diameter                        = getOr(options, 'diameter', 8*options.zoom);
options.clustrules.diameter             = options.diameter;
options.clustrules                      = get_clustrules(options.clustrules);

% Save options back to storage manager
sm.options = options;

% this loads ops1 and checks if processed binary files exist
processed = sm.isRegistrationDone();

% run reg2P if the binaries do not exist
%%%% if tiffs have already been registered, set ops.doRegistration = 0
%%%% and reg2P will just create binary file
% ops1 are the settings and values from registration
if processed
    disp('already registered binary found');
    ops1 = sm.loadRegistrationOptions();
else
    ops1 = reg2P(sm);  % do registration
end

%%
for i = 1:numel(ops1)
    ops         = ops1{i};

    % check if settings are different between ops and options
    % options settings are chosen over ops settings
    ops         = opsChanges(ops, options);

    ops.iplane  = i;

    if numel(ops.yrange)<10 || numel(ops.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
        continue;
    end

    if getOr(ops, {'getSVDcomps'}, 0)
        % extract and write to disk SVD comps (raw data)
        ops    = get_svdcomps(ops, sm);
    end

    if ops.getROIs
        % get sources in stat, and clustering images in res
        [ops, stat, model]           = sourcery(ops, sm);

        if options.fig
            figure(10); clf;
        end

        % extract dF
        switch getOr(ops, 'signalExtraction', 'surround')
            case 'raw'
                [ops, stat, Fcell, FcellNeu] = extractSignalsNoOverlaps(ops, model, stat, sm); %#ok<ASGLU>
            case 'regression'
                [ops, stat, Fcell, FcellNeu] = extractSignals(ops, model, stat, sm); %#ok<ASGLU>
            case 'surround'
                [ops, stat, Fcell, FcellNeu] = extractSignalsSurroundNeuropil2(ops, stat); %#ok<ASGLU>
        end

        % apply user-specific clustrules to infer stat.iscell
        stat                         = classifyROI(stat, ops.clustrules); %#ok<NASGU>


        fileName = sm.getFileForPlane('f', ops.iplane);
        save(fileName,  'ops',  'stat', 'Fcell', 'FcellNeu', '-v7.3');
    end

    fclose('all');
    if ops.DeleteBin
        delete(ops.RegFile);        % delete temporary bin file
    end
end

% clean up
fclose('all');
