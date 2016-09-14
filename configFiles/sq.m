function x = sq(x, varargin)

if numel(varargin)>0
    idims = varargin{1};
    xsiz = size(x);
    idims = idims(xsiz(idims)==1);
    xsiz(idims) = [];
    x = reshape(x, xsiz);
else
    x = squeeze(x);
end

