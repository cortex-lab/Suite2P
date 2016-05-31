function [F0, c, st] = deconvolveL0(F0, kernel, type, Th)
%%
nt0 = numel(kernel);

padKernel = zeros(1,3*nt0-1);
padKernel(nt0 + [1:nt0]) = kernel; 

WtW  = zeros(2*nt0-1, 1);
for k = 1:2*nt0-1
   WtW(k,1) =  padKernel(k+[1:nt0]) * kernel; 
end
%
xrange = [(-nt0+1):(nt0-1)];

F0(end+3*nt0) = 0;
Wf = filter(kernel(end:-1:1), 1, F0);
Wf = Wf(nt0+1:end);

% Wf = conv(F0, kernel);

maxiter = 2000;
c           = zeros(maxiter, 1);
st          = zeros(maxiter, 1);
inner_all   = zeros(maxiter, 1);

NT = numel(Wf);
%%
tots = 0;
icurrent = 0;
while tots<maxiter    
    inner_it = 0;
    while 1 && tots>0     
        % find indices of detected spikes and around them
        inds = repmat(st(1:tots)', 2*nt0-1, 1) + repmat(xrange', 1, tots);
        
        % compute new_coefficients for all new positions
        new_coef = Wf(inds) + WtW * c(1:tots)';
        
        % estimate delta loss upon discarding this spike
        delta = c(1:tots).^2 + 2 * c(1:tots) .* Wf(st(1:tots));
        
        % estimate increase in explained variance from shifting each spike
        switch type
            case 'nonneg'
                DD = max(0, new_coef).^2 - ones(2*nt0-1, 1) * delta';
            case 'real'
                DD = new_coef.^2 - ones(2*nt0-1, 1) * delta';
        end
        
        % find best shift for each spike
        [Mmax, cmax] = max(DD, [], 1);
        
        % find best shift overall
        [M0max, imax] = max(Mmax);
        
        if M0max<3
            % exit the inner loop if it's not worth it
            break;
        else                    
            % add back contribution of shifted spike
            Wf(st(imax) + xrange)  = Wf(st(imax) + xrange)  + c(imax) * WtW;            
            F0(st(imax) + [1:nt0]) = F0(st(imax) + [1:nt0]) + c(imax) * kernel;
            
            % re-assign spike time and magnitude
            c(imax)  = new_coef(cmax(imax), imax);
            st(imax) = st(imax) + cmax(imax) - nt0;           
            
            % remove contribution of shifted spike
            Wf(st(imax) + xrange)  = Wf(st(imax) + xrange) - c(imax) * WtW;
            F0(st(imax) + [1:nt0]) = F0(st(imax) + [1:nt0]) - c(imax) * kernel;
            
            inner_it = inner_it + 1;
        end
    end    
    inner_all(tots+1) = inner_it;    
    
    Cmax = Wf;
    
    Cmax([1:nt0 NT-nt0:NT])= 0;
    
    switch type
        case 'nonneg'
            [Mmax, imax] = max((Cmax(1:NT)));
        case 'real'
            [Mmax, imax] = max(abs(Cmax));
    end
    
    if Mmax<Th
        break;
    end
    
    tots = tots+1;
    icurrent = tots;
    
    c(icurrent) = Wf(imax); 
    st(icurrent) = imax;
    
    
    F0(imax + [1:nt0]) = F0(imax + [1:nt0]) - c(icurrent) * kernel;
    Wf(imax + xrange)        = Wf(imax + xrange) - c(icurrent) * WtW;
    
    
end

c = c(1:icurrent);
st = st(1:icurrent);
