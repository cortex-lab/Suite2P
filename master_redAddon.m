%%
cd('D:\CODE\MariusBox\runSuite2P')
addpath('D:\CODE\GitHub\Suite2P')

% make database to run in batch
make_db_20plane_1x;

ops0.redImageAddon           = 1;

ops0.useGPU                 = 1; % if you can use a GPU in matlab this accelerate registration approx 3 times

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C)
ops0.RegFileRoot            = 'C:/DATA/tempreg';
ops0.RootStorage            = '//zserver4/Data/2P';
ops0.TempStorage            = 'C:/DATA/'; % copy data locally first
ops0.RegFileRoot            = 'C:/DATA/'; 
ops0.ResultsSavePath        = 'D:/DATA/F';

% DO NOT CHANGE: ensures run_pipeline is running correctly
ops0.CopyDataLocally        = 0;
ops0.RegFileTiffLocation    = []; %'D:/DATA/'; % leave empty to NOT save registered tiffs

%%
% ops0.nplanes = 14;
% ops0.nchannels = 2;
ops0.expred = [1 3];
for iexp = 1 % 1:length(db)        %3:length(db)          
   run_REDaddon(iexp, db, ops0) ;
end
%%
