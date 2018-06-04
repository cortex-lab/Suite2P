% rigid registration of frames with offsets ds
function [dreg] = rigidMovie(data, ops1, dsall, yFOVs, xFOVs)

ix0 = 0;
Nbatch = 1000;
dreg = zeros(size(data), 'single');
while ix0<size(data,3)
    indxr = ix0 + (1:Nbatch);
    indxr(indxr>size(data,3)) = [];
    for l = 1:size(xFOVs,2)
        dreg(yFOVs(:,l), xFOVs(:,l), indxr)        = ...
            rigidRegFrames(data(yFOVs(:,l), xFOVs(:,l), indxr), ops1{1,l}, dsall(indxr,:,l));
    end
    ix0 = ix0 + Nbatch;
end