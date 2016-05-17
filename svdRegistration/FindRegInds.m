function ix = FindRegInds(ix,nX,nY)

nA = [nX nY];
ixind = ix > ceil(nA/2);

ix(ixind) = ix(ixind) - (nA(ixind)+1);
ix(~ixind) = ix(~ixind)-1;
