%% set path to your toolbox installation
addpath('C:\CODE\Github\Suite2P\SpikeDetection')
f0      = 30; % sampling frequency
%% load a preoptimized temporal kernel (change this to your toolbox path)
load('D:\CODE\MariusBox\SpikeDetection\kernel.mat')

%% load the procesed data into an array size NT by NN  (timepoints by neurons)

%% pre-process data and extract spikes
% Ffr(econstruction) and Ffp(reprocessed) are outputs

flag_preproc = 1;
[dcell, Ffr, kinterp, Ffp] = run_deconvolution(Ff, f0, kernel, flag_preproc);

%% plots original trace and reconstruction from the model
NN = numel(dcell);
iNN = ceil(rand *NN);
clf
plot(my_conv2(Ff(:,iNN) - dcell{iNN}.B, 2, 1))
hold all
plot(conv(kinterp, Ffr(:,iNN)))
hold off