i = 0;

i = i+1;
db(i).mouse_name    = 'M150116_MD017';
db(i).date          = '2015-02-12';
db(i).expts = [3 5 6 7 9 11 13]; % size tuning
%db(i).expts = [5 6]; % size tuning
db(i).expred = [];
db(i).nchannels = 1;
db(i).gchannel = 1;
db(i).nplanes = 5;
db(i).chunk_align   = 1; %[10 10 10 10 30];
% db(i).planesToProcess = [2:5];

i = i+1;
db(i).mouse_name    = 'M150808_MP016';
db(i).date          = '2015-09-09';
db(i).expts = [5]; % which experiments to cluster together
db(i).expred = [];
db(i).nchannels = 2;
db(i).gchannel = 1;
db(i).nplanes = 5;
db(i).planesToProcess = [2:5];
db(i).regopspath = 'D:\DATA\F\M150808_MP016\2015-09-09\regops_M150808_MP016_2015-09-09_plane';

i = i+1;
db(i).mouse_name    = 'M150127_MD018';
db(i).date          = '2015-02-17';
db(i).expts =[4 5 6 7 9 10 11 12 13]; % which experiments to cluster together
db(i).expred = [];
db(i).nchannels = 2;
db(i).gchannel = 1;
db(i).nplanes = 5;

i = i+1;
db(i).mouse_name    = 'M140924_MD010';
db(i).date          = '2014-10-26';
db(i).expts = [1 3 5 6 8 9 10]; % which experiments to cluster together
db(i).expred = [];
db(i).nchannels = 2;
db(i).gchannel = 1;
db(i).nplanes = 5;


i = i+1;
db(i).mouse_name    = 'M140924_MD010';
db(i).date          = '2014-10-26';
db(i).expts = [3 8 10]; % size tuning
db(i).expred = [];
db(i).nchannels = 2;
db(i).gchannel = 1;
db(i).nplanes = 5;

i = i+1;
db(i).mouse_name    = 'M140909_MD008';
db(i).date          = '2014-09-25';
db(i).expts = [5]; % size tuning
db(i).expred = [];
db(i).nchannels = 2;
db(i).gchannel = 1;
db(i).nplanes = 5;
