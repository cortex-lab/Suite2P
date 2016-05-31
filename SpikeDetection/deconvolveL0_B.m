function [F0, c, st] = deconvolveL0(F0, kernel, type, Th)

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
% Wf(1:nt0) = 0;

% Wf = conv(F0, kernel);

maxiter = 2000;
c           = zeros(maxiter, 1);
st          = zeros(maxiter, 1);
inner_all   = zeros(maxiter, 1);

NT      = numel(Wf);
onesNT1 = ones(NT,1); 
WtOne = sum(kernel)/NT^.5;

C = [1 WtOne; WtOne 1];
Cinv = inv(C);

tic
tots = 0;
while tots<maxiter    
    inner_it = 0;
    while 1 && tots>0        
        inds = repmat(st(1:tots)', 2*nt0-1, 1) + repmat(xrange', 1, tots);
        delta = c(1:tots).^2 + 2 * c(1:tots) .* Wf(st(1:tots));
        new_coef = Wf(inds) + WtW * c(1:tots)';
        switch type
            case 'nonneg'
                DD = max(0, new_coef).^2 - ones(2*nt0-1, 1) * delta';
            case 'real'
                DD = new_coef.^2 - ones(2*nt0-1, 1) * delta';
        end
         
        [Mmax, cmax] = max(DD, [], 1);
        [M0max, imax] = max(Mmax);
        if M0max<10
            break;
        else                    
            Wf(st(imax) + xrange)   = Wf(st(imax) + xrange) + c(imax) * WtW;            
            F0(st(imax) + [1:nt0])  = F0(st(imax) + [1:nt0]) + c(imax) * kernel;
            
            c(imax)             = new_coef(cmax(imax), imax);
            st(imax)            = st(imax) + cmax(imax) - nt0;           
            
            Wf(st(imax) + xrange)  = Wf(st(imax) + xrange)  - c(imax) * WtW;
            F0(st(imax) + [1:nt0]) = F0(st(imax) + [1:nt0]) - c(imax) * kernel;
            
            inner_it = inner_it + 1;
        end
    end    
    inner_all(tots+1) = inner_it;    

    mF0  = sum(F0)/NT^.5;
    coefs = Cinv *      [Wf mF0 * onesNT1]';
    
    Cmax = sum(coefs .* [Wf mF0 * onesNT1]', 1);
    Cmax = Cmax - mF0^2;
    Cmax([1:nt0 NT-nt0:NT])= 0;
    Cmax(coefs(1,:)<0) = 0;
    
    switch type
        case 'nonneg'
            [Mmax, imax] = max(Cmax);
        case 'real'
            [Mmax, imax] = max(Cmax);
    end
    
    if Mmax<Th
        break;
    end
    
    tots = tots+1;
    icurrent = tots;
    
    c(icurrent)     = coefs(1,imax);
    st(icurrent)    = imax;
    
    
    F0(imax + [1:nt0]) = F0(imax + [1:nt0]) - coefs(1, imax) * kernel;
    F0(1:NT) = F0(1:NT) - coefs(2,imax)/NT^.5;
    
    Wf(imax + xrange)  = Wf(imax + xrange)  - coefs(1,imax) * WtW;
    Wf = Wf - coefs(2,imax) * WtOne;
    Wf(1:nt0) = 0;
    
end

c = c(1:icurrent);
st = st(1:icurrent);
