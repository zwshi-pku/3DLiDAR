function show_segs(cloud, labels, gap)
colors = [];
samples = [];
for i = 1:max(labels)
    idxl_sample = (labels == i);
    sample = cloud(idxl_sample,:);
    colors = [colors; repmat(rand(1,3),size(sample,1),1)];
    samples = [samples; sample];
end
pcshow(samples(1:gap:end,:),colors(1:gap:end,:))
end