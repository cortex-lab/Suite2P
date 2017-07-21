% shifts indices to +/- rather than circular coordinates
%input   1 2 ... nA-1 nA
% ->     0 1 ...  -2  -1
function ix = FindRegInds(ix,nA)

nA    = repmat(nA, size(ix,1), 1);

ixind = ix > ceil(nA/2);

ix(ixind) = ix(ixind) - (nA(ixind)+1);
ix(~ixind) = ix(~ixind)-1;
