function Smooth = expfilt(S1, sig)

NN = size(S1,1);
NT = size(S1,2);

dt = 1:1:(5*sig);
filt = exp( - dt/sig);
filt = filt'/sum(filt);

% Norms = conv(ones(NT,1), gaus, 'same');
%Smooth = zeros(NN, NT);
%for n = 1:NN
%    Smooth(n,:) = (conv(S1(n,:)', gaus, 'same')./Norms)';
%end

Smooth = filter(filt, 1, [zeros(4*sig, NN+1);S1' ones(NT,1)]);
Smooth = Smooth(1+4*sig:end, :);
Smooth = Smooth(:,1:NN) ./ (Smooth(:, NN+1) * ones(1,NN));

Smooth = Smooth';