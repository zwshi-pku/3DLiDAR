%% Fusion power line
function [powerLines_pro, index]= merge(powerLines, ind)
for i=1:size(ind,1)
    if powerLines(ind(i)).Label == 1
        continue;
    end
    ids = findMerge(powerLines, ind, i);
    powerLines(ind(i)).Ids = ids;
    for j=1:size(ids,1)
        powerLines(ids(j)).Label = 1;        
    end
end

% Delete marked data
powerLines_pro = powerLines;
i = 1;
while i <= size(powerLines_pro,2)
    if powerLines_pro(i).Label
        powerLines_pro(i) = [];
    else
        i = i + 1;
    end
end

for i = 1:size(powerLines_pro,2)
    for j =1:size(powerLines_pro(i).Ids,1)
        powerLines_pro(i).Location = [powerLines_pro(i).Location; powerLines(powerLines_pro(i).Ids(j)).Location];
    end    
end

counts = [];
for i = 1:size(powerLines_pro,2)
    powerLines_pro(i).Label = 0;
    powerLines_pro(i).Count = size(powerLines_pro(i).Location,1);
    powerLines_pro(i).Ids = [];
    counts = [counts; powerLines(i).Count];
end
[counts_new, index] = sort(counts,'descend');

end