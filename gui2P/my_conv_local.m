
function Smooth = my_conv_local(S1, sig)
NN = size(S1,1);
NT = size(S1,2);
dt = -4*sig:1:4*sig;
gaus = exp( - dt.^2/(2*sig^2));
gaus = gaus'/sum(gaus);
Smooth = filter(gaus, 1, [S1' ones(NT,1); zeros(4*sig, NN+1)]);
Smooth = Smooth(1+4*sig:end, :);
Smooth = Smooth(:,1:NN) ./ (Smooth(:, NN+1) * ones(1,NN));
Smooth = Smooth';