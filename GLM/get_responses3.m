function [R, filt, STA, frec] = get_responses3(F0, ton, ops, varargin)

filt = [];
frec = [];
switch ops.method
    case 'deconv'
        if ~isempty(varargin)
            [R, filt, frec] = deconv(F0, ton, ops,varargin{1});
        else
            [R, filt, frec] = deconv(F0, ton, ops);
        end
        [~, ~, STA] = get_sta(F0, ton, ops);

    case 'sta'
        if ~isempty(varargin)
            [R, filt, STA] = get_sta(F0, ton, ops, varargin{1});
        else
             [R, filt, STA] = get_sta(F0, ton, ops);
        end
end

end

function [R, filt, frec] = deconv(F0, ton, ops, varargin)


% initialize filters
if numel(varargin)<1
    filt = [0 0 exp(-[1:1:(ops.ntf-2)]/(ops.ntf/3))];
    filt = filt/sum(abs(filt(:)));
else
    filt = varargin{1};
    ops.niter = 1;
end
[NN, NT] = size(F0);

for i = 1:length(ton)
   ton{i}(ton{i}+ops.ntf>NT) = []; 
end


Nstims = numel(ton);
stims = zeros(Nstims,NT);
for i = 1:Nstims
    stims(i,ton{i}) = 1;
end
if ops.smoothORI
   stims(:, :)= my_conv_circ(stims(:, :)', ops.smoothORI)';    
end

pred2 = zeros(ops.ntf, NT);

F1 = sum(F0,1);

for n = 1:ops.niter
    
    pred = filter(filt', 1, stims')';
        
    X = [pred; ones(1,NT)];
    xtx = (X*X')/NT;
    f0x = F0*X'/NT;
    Lreg = eye(size(xtx,1));
    Lreg(Nstims+1:end, Nstims+1:end) = 0;
    B = f0x/(xtx + ops.lam * Lreg);
    %      vexp = mean(mean((F0 - B *X).^2)/VF;
    %     fprintf('error is %2.5f\n', vexp)
    
    if n ==ops.niter
       frec = B *X; 
    end
    if numel(varargin)<1
        for i = 1:ops.ntf
            for j = 1:Nstims
                pred2(i, ton{j}+i-1) = sum(B(:,j));
            end
        end
        pred2(:, NT+1:end) = [];
        
        X = [pred2; ones(1,NT);];
        xtx = (X*X')/NT;
        f0x = F1*X'/NT;
        B2 = f0x/xtx;
        filt = B2(1:ops.ntf);
        filt = filt/sum(abs(filt));
    end
       
%     fprintf('error is %2.5f\n', vexp)
end


R = B(1:NN,1:Nstims);

end

function [R, U, Fresp] = get_sta(F0, ton, ops, varargin)
dt = -10:1:60;
[Npops NT] = size(F0);

Nstims = numel(ton);

Fresp = zeros(Npops, Nstims, numel(dt));
for i = 1:Nstims
    tclipped = ton{i};
    tclipped(tclipped+dt(end)>NT | tclipped+dt(1)<1) = [];
    ts = repmat(tclipped, numel(dt), 1) + repmat(dt', 1, numel(tclipped));
   
    Fresp(:,i, :) = mean(reshape(F0(:,ts), [size(F0,1) size(ts)]), 3);
end


%%
% Fbase = mean(mean(Fresp(:,:,1:60),3),2);
% Fresp = Fresp - repmat(Fbase, 1, size(Fresp,2), size(Fresp,3));

Fresp = reshape(Fresp, [], numel(dt));
if ~isempty(varargin)
   U = varargin{1};
else
    [U Sv] = eigs(double(Fresp' * Fresp), 1);
    [~, imax] = max(abs(U));
    U = U * sign(U(imax));
end

% R = mean(Fresp, 2);
R = Fresp * U(:);

R = reshape(R, Npops, Nstims);

Fresp = reshape(Fresp, Npops, Nstims, numel(dt));
end
