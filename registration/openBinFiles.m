% open binary files to write registered tiffs to
% takes ops, ops1, and red_binary (if 1 then write red channel to binary)
function [ops1, fid, fidRED, fidIntpol] = openBinFiles(sm, ops1, red_binary)
numPlanes = numel(sm.options.planesToProcess);
fid    = cell(numPlanes, size(ops1, 2));
fidRED = cell(numPlanes, size(ops1, 2));
fidIntpol = [];
regFileBinLocation = sm.getOr('RegFileBinLocation');
doInterpolation = sm.getOr('interpolateAcrossPlanes') && ~isempty(regFileBinLocation);
if doInterpolation
    fidIntpol = cell(numPlanes, size(ops1,2));
end

for i = 1:numPlanes
    for j = 1:size(ops1, 2)
        planeInd = i + (j-1)*numPlanes;
        ops1{i, j}.RegFile = sm.registrationFileOnFastStorage(planeInd);

        % open bin file for writing
        fid{i, j} = fopen(ops1{i, j}.RegFile, 'w');

        if red_binary
            ops1{i, j}.RegFile2 = sm.registrationFileOnFastStorage(planeInd, 'red');
            fidRED{i, j} = fopen(ops1{i, j}.RegFile2, 'w');
        end

        if doInterpolation
            % open separate files for result after averaging across
            % neighbouring planes
            fidIntpol{i, j} = fopen(sm.registrationFileOnSlowStorage(planeInd, 'interp'), 'w');
        end
    end
end