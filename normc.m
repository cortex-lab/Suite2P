function v = normc(v)

v = v./repmat(sum(v.^2, 1), [size(v,1),ones(1,ndims(v)-1)]).^.5;