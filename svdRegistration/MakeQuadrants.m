function [ops,numBlocks] = MakeQuadrants(ops,yB,xB)
 
ib = 0;
for iy = 1:length(yB)-1
  for ix = 1:length(xB)-1
    ib = ib+1;
    ops.yBL{ib} = (yB(iy)+1):yB(iy+1);
    ops.xBL{ib} = (xB(ix)+1):xB(ix+1);
  end
end
numBlocks = ib;
