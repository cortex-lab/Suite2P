function B = nDshift(B, sig, varargin)

idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end
if numel(idims)>1 && numel(sig)>1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

for i = 1:length(idims)
    sig = sigall(i);
    
    idim = idims(i);
    Nd = ndims(B);
    
    B = permute(B, [idim 1:idim-1 idim+1:Nd]);
    
    dsnew = size(B);
    
    B = reshape(B, dsnew(1), []);
    dsnew2 = size(B);
    
    if sig>0
        irange = sig+1:dsnew2(1);
        B(irange, :) = B(irange-sig, :);
        B(1:sig, :) = 0;
    else
        sig = -sig;
        irange = 1:dsnew2(1)-sig;
        B(irange, :) = B(irange+sig, :);
        B([1:sig] + dsnew2(1)-sig, :) = 0;
    end    
    
     B = reshape(B, dsnew);
       
    B = permute(B, [2:idim 1 idim+1:Nd]);
end