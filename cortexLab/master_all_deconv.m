addpath('D:\CODE\MariusBox\SpikeDetection')

f0      = 30; % sampling frequencys

mname   = 'M150329_MP009';
datExp  = '2015-04-27';
root    = 'D:\DATA\F\';
plane_um = 0;
pix_um  = 2;
npl     = 1;

%% load the procesed data
[Ff, res] = loadPROC(root, mname, datExp, npl, plane_um, pix_um, f0);


%% ITeRATIVE OPTIMIZATION OF PARAMETERS
% load('D:\CODE\MariusBox\BigNeuralCode\results\driv1112F.mat')
load('D:\CODE\MariusBox\SpikeDetection\kernel.mat')

[dcell, Ffr, kinterp] = run_deconvolution(Ff, f0, kernel);

%%
NN = numel(dcell);
iNN = ceil(rand *NN);
clf
plot(my_conv2(Ff(:,iNN) - dcell{iNN}.B, 2, 1))
hold all
plot(conv(kinterp, Ffr(:,iNN)))
hold off