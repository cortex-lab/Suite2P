function run_REDaddon_sourcery(db, ops0)
ops = build_ops3(db, ops0);
% red channel addon to already processed data
mimgG = [];

if sum(ismember(ops.expts, ops.expred))>0
    mimgR = regRedChannelExpts(ops);
else
    if ops.nchannels_red==2
        % return green mean image from short red/green recording
        [mimgR, mimgG] = regRedGreenChannel(ops);
    else
        mimgR = regRedChannelOnly(ops);
    end
end

save(fullfile(ops.ResultsSavePath,'redchannel.mat'), 'mimgR','mimgG');

add_red_channel_sourcery(ops);
