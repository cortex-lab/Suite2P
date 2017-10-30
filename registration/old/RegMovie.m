% registers frames using offsets dsall (rigid registration)
% loops over batches of frames and over splits in frames
function [dreg] = RegMovie(data, ops1, dsall, yFOVs, xFOVs)

dreg = zeros(size(data), 'like', data);
for l = 1:size(xFOVs,2)
  dreg(yFOVs(:,l), xFOVs(:,l), indxr)        = ...
      rigidRegMovie(data(yFOVs(:,l), xFOVs(:,l), :), ops1{1,l}, dsall(:,:,l));
end
  