function ops = buildRegOps(ops)

% number of channels in db.expts that aren't db.expred
ops.nchannels      = getOr(ops, 'nchannels', 1);
ops.nplanes        = getOr(ops, 'nplanes', 1);

% --- registration options --- %
ops.kriging         = getOr(ops, 'kriging', 1);  % kriging subpixel alignment
ops.PhaseCorrelation = getOr(ops, 'PhaseCorrelation', 1); % set to 0 for non-whitened cross-correlation
if ~ops.kriging
    ops.registrationUpsample = getOr(ops, {'registrationUpsample'}, 1);  % upsampling factor during registration
    % 1 for no upsampling is much faster, 2 may give better subpixel accuracy
end
ops.NimgFirstRegistration  = getOr(ops, 'NimgFirstRegistration', 500); % number of images to include in the first registration pass
ops.NiterPrealign   = getOr(ops, {'NiterPrealign'}, 10); % number of iterations for initial mean image
ops.showTargetRegistration = getOr(ops, 'showTargetRegistration', 1); % shows the image targets for all planes to be registered

% --- non-rigid registration splits FOV into multiple blocks, registers
% each block separately, and then smooths shifts over blocks --- %
% numBlocks = [y-blocks x-blocks]
if getOr(ops, 'nonrigid', 0)
    % default is 8 blocks in the Y direction (1/6 pixels each)
    % ops.numBlocks(1) = # of blocks in Y, ops.numBlocks(2) = # of blocks in X
    ops.numBlocks      = getOr(ops, {'numBlocks'}, [8 1]);
    disp('non-rigid registration chosen');
elseif isfield(ops, 'numBlocks') && ~isempty(ops.numBlocks)
    % if numBlocks is set to a value greater than 1, then run non-rigid
    if numel(ops.numBlocks) == 1
        ops.numBlocks = [ops.numBlocks 1];
    end
    if sum(ops.numBlocks) > 2
        ops.nonrigid               = 1;
        disp('non-rigid registration chosen');
    end
else
    ops.numBlocks = [1 1];
    ops.nonrigid = 0;
end

% --- keep frames at beginning and end of expts (to check for movement)
ops.nimgbegend     = getOr(ops, {'nimgbegend'}, 0);

% --- write tiffs to RegFileTiffLocation if field is not empty --- %
ops.RegFileTiffLocation = getOr(ops, {'RegFileTiffLocation'}, []);
% --- copy binary file to regfilebinlocation if field is not empty --- %
ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);

% --- for aligning across planes --- %
ops.alignTargetImages  = getOr(ops, {'alignTargetImages'}, false); % if true, align target images to each other
ops.interpolateAcrossPlanes = getOr(ops, {'interpolateAcrossPlanes'}, false); %if true, similar looking planes will be averaged together to generate final movie
ops.planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:ops.nplanes); % these planes will be considered for interpolation
ops.alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false); % at each time point, frame will be aligned to best matching target image (from different planes)

% --- for within plane splitting (if tiffs are too large in Y and X) --- %
ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]); % split FOV into chunks if memory issue
% ops.splitFOV(1) = # of subsets in Y, ops.splitFOV(2) = # of subsets in X
ops.smooth_time_space  = getOr(ops, 'smooth_time_space', []);

% --- bidirectional phase offset computation --- %
% assumes same bidirectional offset for each plane
ops.dobidi             = getOr(ops, {'dobidi'}, 1); % compute bidiphase?
% if set to a value by user, do not recompute
if isfield(ops, 'BiDiPhase')
    ops.dobidi         = 0;
end
ops.BiDiPhase          = getOr(ops, {'BiDiPhase'}, 0); % set to default 0