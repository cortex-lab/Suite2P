function ops = addfields(ops, opsnew)

fieldNames = fieldnames(opsnew);
for j = 1:size(fieldNames,1)
    if ~isempty(opsnew.(fieldNames{j}))
        ops.(fieldNames{j}) = opsnew.(fieldNames{j});
    end
end