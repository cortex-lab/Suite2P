function [B err] = get_neurop(Ff, F1, kerneli, npad, Fneu)

X = conv(kerneli, F1);

X = X([1:numel(F1)]);

if nargin<5
    A = cat(2, X, ones(numel(F1), 1));
else
    A = cat(2, X, ones(numel(F1), 1), Fneu);
end

AtA = (A'*A) / size(A,1) + 1e-6;
FtA = (Ff'*A) / size(A,1);

B = FtA / AtA;

err = std(Ff - A * B');