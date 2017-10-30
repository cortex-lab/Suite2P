% function to overwrite ops1 fields with ops0 fields if settings have
% changed for new runs
function ops = opsChanges(ops1, ops0)

ops = ops1;
f = fieldnames(ops0);

for j = 1:length(f)
    ops.(f{j}) = ops0.(f{j});
end
