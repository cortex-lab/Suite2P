load examp_dset

%% high-pass filter the data
ops.filtering_method = 'median'; % options are median, exp or none

switch ops.filtering_method
    case 'median'
        ops.medwind         = 60; % window of filtering in timepoints
        Fm = F;
        for i = 1:size(F,1)
            Fm(i,:) = F(i,:) - fastmedfilt1d(F(i,:)', ops.medwind)';
        end
    case 'exp'
        ops.expwind = 30; % mean of exponential filter in timepoints
        Fm = F - expfilt(F, ops.expwind);
    case 'none'
        Fm  = F;        
end

%% separate stimuli into two sets for crossvalidation
for i = 1:length(ton)
    r = randperm(numel(ton{i}));
    ton1{i} = ton{i}(r(1:ceil(numel(r)/2)));
    ton2{i} = ton{i}(r(1+ceil(numel(r)/2):numel(r)));
end

%%
ops.ntf         = 60; % timepoints for the filter
ops.niter       = 10; % number of iterations to refine the filter

ops.method      = 'deconv'; % deconv or sta (stimulus-triggered average)

ops.lam             = 1e-5; % regularizer for regression problem
ops.smoothORI       = 0;    % smoothing width (in orientation steps), 0 is no smoothing

[R1, filt1]           = get_responses3(Fm, ton1, ops, filt);
[R2, filt2]           = get_responses3(Fm, ton2, ops, filt);

[bestphi1, amps1, widths1] = get_vonMisesFit(R1);
[bestphi2, amps2, widths2] = get_vonMisesFit(R2);

NN              = numel(bestphi1);
delta_phi       = rem(abs(bestphi1 - bestphi2), pi)*180/pi;
delta_phi_shuff = rem(abs(bestphi1 - bestphi1(randperm(NN))), pi)*180/pi;

fprintf('Mean orientation error %2.2f degrees\n', mean(delta_phi));
fprintf('Mean shuffled orientation error %2.2f degrees\n', mean(delta_phi_shuff));

%%
