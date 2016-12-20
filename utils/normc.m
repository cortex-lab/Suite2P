function v = normc(v, eps, ndim)
% add epsilon to normalizer, and choose what dimension ndim to normalize

if nargin>2
    v = permute(v, [ndim 1:ndim-1 ndim+1:ndims(v)]);
    dimv = size(v);
    v = v(:,:);
end

norms = sum(v.^2, 1).^.5;
if nargin>1
    norms = norms + eps;
end

if nargin==1
    v = bsxfun(@rdivide, v, norms);
end
if nargin>2
    v = reshape(v, dimv);
    v = permute(v, [2:ndim 1 ndim+1:ndims(v)]);
end

