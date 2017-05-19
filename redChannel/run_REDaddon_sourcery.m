function run_REDaddon_sourcery(iexp, db, ops0)
% red channel addon to already processed data
mimgG = [];
ops = build_ops3(db(iexp), ops0);

if sum(ops.expts==ops.expred)>0
    mimgR = red_channel_mean(ops);
else
    if ops.nchannels_red==2
        % return green mean image from recording
        [mimgR,mimgG] = red_channel_mean3(ops);
    else
        mimgR = red_channel_mean2(ops);
    end
end

save(fullfile(ops.ResultsSavePath,'redchannel.mat'), 'mimgR','mimgG');

add_red_channel_sourcery(ops);
