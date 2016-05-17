function icl = get_cl(iclust, Nk)

icl = cell(Nk,1);

for i = 1:Nk
    ix = find(iclust==i);
    icl{i} = ix;
end
