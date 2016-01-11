function run_REDaddon(iexp, db, ops0)
% red channel addon to already processed data

ops = build_ops(db(iexp), ops0);
mimgR = red_channel_mean(ops);
add_red_channel(ops, mimgR);
