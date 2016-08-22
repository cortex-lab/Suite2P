function dSplit = divide_dcell(dcell, sz)
%%
dSplit = cell(numel(sz), 1);

offset = 0;
for i = 1:numel(sz);
    last_point = offset + sz(i);
    dSplit{i}.dcell = dcell;
    for j = 1:numel(dcell)
        its = dcell{j}.st>offset & dcell{j}.st<=last_point;
        dSplit{j}.dcell{j}.c = dcell{j}.c(its);
        dSplit{j}.dcell{j}.st = dcell{j}.st(its) - offset;
    end
    
    
    offset = last_point;
end