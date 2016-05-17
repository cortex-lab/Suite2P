function ix = InvRegInds(ix,nX,nY)

nA = [nX nY];
ixind = ix < 0;

ix(ixind) = ix(ixind) + (nA(ixind)+1);
ix(~ixind) = ix(~ixind) + 1;

