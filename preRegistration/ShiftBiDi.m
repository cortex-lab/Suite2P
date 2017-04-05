function IMG = ShiftBiDi(BiDiPhase, IMG, Ly, Lx)

yrange = 2:2:Ly;
if BiDiPhase>0
    IMG(yrange,(1+BiDiPhase):Lx,:,:) = IMG(yrange, 1:(Lx-BiDiPhase),:,:);
else
    IMG(yrange,1:Lx+BiDiPhase,:,:)   = IMG(yrange, 1-BiDiPhase:Lx,:,:);
end