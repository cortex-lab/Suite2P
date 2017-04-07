function [w_out inomax]= my_ica(ca, neu, fs, w_default, w_max)
%%
if nargin<4; w_default = .7; end
if nargin<5; w_max = 1.4; end
    
ca = double(ca);
neu = double(neu);

[b1, a1] = butter(3, 2*[1/100 1/10]/fs, 'bandpass');

nPerBatch = 5e8/size(ca,1)/8;

for j = 1:ceil(size(ca,2)/nPerBatch)
    i0 = (j-1) * nPerBatch + [1:nPerBatch];
    i0(i0>size(ca,2)) = [];
    
    ca(:,i0)   = filtfilt(b1, a1, ca(:,i0));
    neu(:,i0) = filtfilt(b1, a1, neu(:,i0));
    
end

neu = bsxfun(@minus, neu, mean(neu,1));
ca   = bsxfun(@minus, ca, mean(ca,1));
    
wR = mean(ca.*neu, 1)./mean(neu.*neu, 1);
ca = ca - bsxfun(@times, neu, wR);

[NT NN] = size(ca);
ws = linspace(-1, 1, 21);

sk = zeros(size(ca,2), length(ws));
for i = 1:length(ws)
    dF = ca - ws(i) * neu;
    sk(:,i) = skewness(dF, 1, 1);
end

% find all the local maxima
pks = sk(:, 2:end-1)> sk(:, 1:end-2) & sk(:, 2:end-1)> sk(:, 3:end);
pks = bsxfun(@rdivide, pks , (1 + abs(ws(2:end-1))));
[Mmax, imax] = max(pks, [], 2);
inomax = Mmax<1e-10; 

w1 = ws(1+imax);

% refine the grid
ws = linspace(-.1, .1, 21);
for i = 1:length(ws)
    c = bsxfun(@plus, w1, ws(i));
    dF = ca - bsxfun(@times, neu, c);
    sk(:,i) = skewness(dF, 1, 1);
end
[Mmax, imax] = max(sk, [], 2);
w2 = w1 + ws(imax);

w_out = w2 + wR;
w_out(inomax) = w_default;

w_out = w_out(:);
w_out = min(w_max, w_out);
% for iter = 1:niter
%    % compute projection
%    
%    wtz = bsxfun(@times, ca, w(1,:)) + bsxfun(@times, neu, w(2,:));
%    
%    w(1,:) = mean(ca  .* wtz.^3 - 3*bsxfun(@times, wtz.^2, w(1,:)), 1);
%    w(2,:) = mean(neu .* wtz.^3 - 3*bsxfun(@times, wtz.^2, w(2,:)), 1);
%    
%    w = normc(w);
%    w = bsxfun(@times, w, sign(w(1,:)));
% end

% w_out = -w(2,:)./w(1,:) + wR;
% 


