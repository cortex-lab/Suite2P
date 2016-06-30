function run_REDaddon(iexp, db, ops0)
% red channel addon to already processed data

ops = build_ops3(db(iexp), ops0);

if sum(ops.expts==ops.expred)>0
    mimgR = red_channel_mean(ops);
else
    if ops.nchannels_red==2
        mimgR = red_channel_mean3(ops);
    else
        mimgR = red_channel_mean2(ops);
    end
end

add_red_channel(ops, mimgR);
