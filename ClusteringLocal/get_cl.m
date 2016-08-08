function icl         = get_cl(iclust, Nk)

icl = cell(Nk,1);

for k = 1:Nk
  icl{k} = find(iclust==k);  
end
