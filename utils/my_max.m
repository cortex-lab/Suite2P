function S1 = my_max(S1, sig, varargin)

if numel(varargin)>0
    S1 = -my_min(-S1, sig, varargin{1});
else
    S1 = -my_min(-S1, sig);
end