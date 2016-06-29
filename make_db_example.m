i = 0;

i = i+1;
db(i).mouse_name    = 'M150329_MP009';
db(i).date          = '2015-04-29';
db(i).expts         = [4 5 6];
%
i = i+1;
db(i).mouse_name    = 'M150329_MP009';
db(i).date          = '2015-04-10';
db(i).expts         = [5 6 7 8 9 10 11];

i = i+1;
db(i).mouse_name    = 'M150824_MP019';
db(i).date          = '2015-12-19';
db(i).expts         = [4];

% example extra entries
% db(i).AlignToRedChannel= 1;
% db(i).BiDiPhase        = 0; % adjust the relative phase of consecutive lines
% db(i).nSVD             = 1000; % will overwrite the default, only for this dataset
% db(i).comments      = 'this was an adaptation experiment';
% db(i).expred        = [4]; % say one block which had a red channel 
% db(i).nchannels_red = 2; % how many channels did the red block have in total (assumes red is last)