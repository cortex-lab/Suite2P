
function I = redraw_figure(h)
% 
Sat1     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
Sat2     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
H1              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
H2              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);

[iclust1, iclust2, V1, V2] = ...
    getviclust(h.dat.stat, h.dat.cl.Ly,  h.dat.cl.Lx, h.dat.cl.vmap, h.dat.F.ichosen);

iselect     = iclust1==h.dat.F.ichosen;
Sat1(iselect)= 0;

iselect     = iclust2==h.dat.F.ichosen;
Sat2(iselect)= 0;

H1(iclust1>0)   = h.dat.cl.rands(iclust1(iclust1>0));
H2(iclust2>0)   = h.dat.cl.rands(iclust2(iclust2>0));

I = hsv2rgb(cat(3, H1, Sat1, V1));
I = min(I, 1);
axes(h.axes2); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off

I = hsv2rgb(cat(3, H2, Sat2, V2));
I = min(I, 1);
axes(h.axes3); imagesc(I);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off

end

function [iclust1, iclust2, V1, V2] = getviclust(stat, Ly, Lx, vmap, ichosen)

iclust1 = zeros(Ly, Lx);
iclust2 = zeros(Ly, Lx);
V1      = zeros(Ly, Lx);
V2      = zeros(Ly, Lx);

for j = 1:numel(stat)
    ipix    = stat(j).ipix;
    lambda   = stat(j).lambda;
    
    if ichosen==j
        inew = true(numel(ipix), 1);
    else
        if stat(j).iscell
            inew    = lambda(:)>(V1(ipix) + 1e-6);
        else
            inew    = lambda(:)>(V2(ipix) + 1e-6);
        end
    end
    
    switch vmap
        case 'var'
            L0      = stat(j).lambda(inew);
        case 'unit'
            L0      = stat(j).lam(inew);    
    end
    if stat(j).iscell
        V1(ipix(inew))      = L0;
        iclust1(ipix(inew)) = j;
    else
        V2(ipix(inew))      = L0;
        iclust2(ipix(inew)) = j;
    end
end
mV = mean([V1(V1>0); V2(V2>0)]);
V1 = V1/mV;
V2 = V2/mV;
end