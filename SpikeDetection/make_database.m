
cd('D:\DATA\Chen203\processed_data')
files = dir('data*.mat');

clear Nspikes
for i = 1:length(files)
    load(files(i).name)
    % load('data_20120521_cell4_003.mat');
    %---------------------------------------

    dat(i).fmean_roi=obj.timeSeriesArrayHash.value{1}.valueMatrix;
    dat(i).fmean_neuropil=obj.timeSeriesArrayHash.value{2}.valueMatrix;
%     fmean_comp=fmean_roi-0.7*fmean_neuropil;
    dat(i).t_frame=obj.timeSeriesArrayHash.value{1}.time;
    
    filt=obj.timeSeriesArrayHash.value{4}.valueMatrix;
    t_ephys=obj.timeSeriesArrayHash.value{4}.time;
    
    detected_spikes=obj.timeSeriesArrayHash.value{5}.valueMatrix;
    dat(i).spike_time=t_ephys(detected_spikes);
    
    Nspikes(i) = sum(detected_spikes);
    switch obj.indicator
        case 'GCaMP6f'
            dat(i).type = 1;
        case 'GCaMP6s'
            dat(i).type = 2;
        case 'GCaMP5'
            dat(i).type = 3;
    end
    dat(i).name = files(i).name;
end

%%
for i = 1:length(files)
    Fneu = dat(i).fmean_roi - .7 * dat(i).fmean_neuropil;
    NT = size(Fneu,1);
    
    Fneu = cat(1, Fneu(100:-1:1), Fneu, Fneu(end:-1:end-100));
    
    Fneu = Fneu - ordfilt2(Fneu, 31, true(201,1));
    Fneu = Fneu(100 + [1:NT]);
    
%     Fsort = sort(Fneu, 'ascend');
%     Fneu = Fneu - Fsort(ceil(1/10 * numel(Fneu)));
        
    dat(i).Fneu = Fneu;
    
%     plot(Fneu)
%     drawnow 
    
    if i==12
        dat(i).quality =0;
    else
         dat(i).quality =1;
    end
end
%% put into single cells
dcell = [];
clear fname_list
for i = 1:length(dat)
    fname_list{i} = dat(i).name(1:end-4-3);
end
[~, ~, IC] = unique(fname_list);
%
for icell = 1:max(IC)
    dcell{icell}.F = [];
    dcell{icell}.tspikes = [];
    totframes = 0;
    ix = 1;
    tfr{icell} = 0;
    for i = 1:length(dat)
        if dat(i).quality && IC(i)==icell
            Fneu = dat(i).Fneu;
           dcell{icell}.F = cat(1, dcell{icell}.F, Fneu);
           [~, idx] = histc(dat(i).spike_time, dat(i).t_frame);
           dcell{icell}.tspikes = cat(1, dcell{icell}.tspikes, idx + totframes);
           totframes = totframes + numel(Fneu);
           
           ix = ix + 1;
           tfr{icell}(ix) = totframes;
        end
    end
    
    maxtsp = max(dcell{icell}.tspikes);
    dcell{icell}.F(maxtsp+220:end) = [];
end







