function run_REDaddon(iexp, db, ops0)
% red channel addon to already processed data
mimgG = [];
ops = build_ops3(db(iexp), ops0);

if sum(ismember(ops.expts, ops.expred)) == numel(ops.expts)
    mimgR = red_channel_mean(ops);
else
    if ops.nchannels_red==2
        % return green mean image from recording
        [mimgR,mimgG] = red_channel_mean3(ops);
    else
        mimgR = red_channel_mean2(ops);
    end
end

add_red_channel(ops, mimgR, mimgG);
