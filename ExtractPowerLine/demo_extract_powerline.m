% /****************************************************************************
%  * Please cite the following paper, If you use this code in your work.
%  *
%  * Zhenwei Shi, Yi Lin, and Hui Li. "Extraction of urban power lines and 
%  * potential hazard analysis from mobile laser scanning point clouds." 
%  * International Journal of Remote Sensing 41, no. 9 (2020): 3411-3428.
%  *
%  * The paper can be downloaded from
%  * https://www.tandfonline.com/doi/full/10.1080/01431161.2019.1701726
%  * Copyright
%  * Institute of Remote Sensing & GIS, 
%  * School of Earth and Space Sciences, Peking University (www.pku.edu.cn)
%  * Zhenwei Shi; Yi Lin; Hui Li
%  * contact us: zwshi@pku.edu.cn; lihui@pku.edu.cn
% *****************************************************************************/

% Extracting powerline from mobile lidar point cloud
% function: [isPLIndex] = extractPLs(pointcloud,radius,angleThr,LThr)
% pointcloud: mobile lidar point cloud [x y z]
% radius: Neighborhood point search radius
% angleThr: Threshold value of the angle between normal vector and (0 0 1)
% LThr: linear feature L satisfy the conditions L≥LThr for power line
% point.

% Example
% radius = 0.5;
% angleThr = 10;
% LThr = 0.98;
% [isPLIndex] = extractPLs(pointcloud,radius,angleThr,LThr);

%% step1 Mobile LiDAR filtering
clc; clear; close all;
tic
nonGroundPoints = [];
for i=8:10
    filename = ['L037_Sens1_600x250_cloud_' num2str(i)];
    load([filename '.mat'])
    
    [counts,centers] = hist(ptCloudA.data(:,3),50);
    [~,I] = max(counts);
    nonGroundIndex = ptCloudA.data(:,3)>centers(I+3);  
    nonGroundPoints = [nonGroundPoints; ptCloudA.data(nonGroundIndex,1:3)];
end
toc
figure
pcshow(nonGroundPoints)
title('Non ground points')
view(60,25)
print(gcf,'-dpng','-r300', 'f1_nonGroundPoints.png')

%% step1 Extract powerline point points
radius = 0.5;
angleThr = 10;
LThr = 0.98;
tic
isPLIndex = extractPLs(nonGroundPoints,radius,angleThr,LThr);
toc
isPLIndex = logical(isPLIndex);

figure
pcshow(nonGroundPoints(isPLIndex,:))
title('Candidate power line points')
view(60,25)
print(gcf,'-dpng','-r300', 'f2_candidate powerline points.png')

%% step2 Euclidean clustering
ptCloud = pointCloud(nonGroundPoints(isPLIndex,:));
minDistance = 2.0;
[labels,numClusters] = pcsegdist(ptCloud,minDistance);

figure
pcshow(ptCloud.Location,labels)
title('Candidate power line clusters')
view(60,25)
print(gcf,'-dpng','-r300', 'f3_candidate powerline points clusters.png')

index = ones(size(labels,1),1);
power_lines = [];
for i=1:numClusters
    cluster = index*i == labels;
    if(sum(cluster)<15)
       continue; 
    end    
    xyzs = ptCloud.Location(cluster,:);
    power_lines = [power_lines; xyzs];
end
PLs = pointCloud(power_lines);
figure
pcshow(PLs);
title('Power line clusters')
view(60,25)
print(gcf,'-dpng','-r300', 'f4_powerline points clusters.png')

%% step3 Power line modeling
ptCloud = PLs;
minDistance = 0.3;
[labels,numClusters] = pcsegdist(ptCloud,minDistance);
figure
show_segs(ptCloud.Location,labels,1);
title('Different clusters are given different colors')
view(60,25)
print(gcf,'-dpng','-r300', 'f5_colorization clusters.png')

counts = [];
index = ones(size(labels,1),1);

% Constructing structure based on Clustering
for i=1:numClusters
    cluster_index = index*i == labels;
    cluster_raw = ptCloud.Location(cluster_index,:); 
    xyzs_new = cluster_raw;
    powerLines(i).Location = xyzs_new;
    powerLines(i).Label = 0;
    powerLines(i).Count = size(xyzs_new,1);
    powerLines(i).Ids = [];
    counts = [counts; powerLines(i).Count];
end

% Get the sorting ID of cluster points
[counts_new, ind] = sort(counts,'descend');
[powerLines_pro ind]= merge(powerLines, ind);
[powerLines_pro ind]= merge(powerLines_pro, ind);
[powerLines_pro ind]= merge(powerLines_pro, ind);

powerLines_new = [];
colors = [];
for i = 1:size(powerLines_pro,2)
    powerLines_new = [powerLines_new; powerLines_pro(i).Location];
    colors = [colors; repmat(rand(1,3),size(powerLines_pro(i).Location,1),1)];
end
% ptpl = pointCloud(powerLines_new,'Color',colors)
figure
% pcshow(ptpl)
pcshow(powerLines_new,colors)
title('Power line clusters')
view(60,25)
print(gcf,'-dpng','-r300', 'f6_powerLines clusters.png')

% Calculate the power line length and delete the power line less than the length threshold
i = 1;
while i <= size(powerLines_pro,2)
    dist = getDist(powerLines_pro(i).Location);
    if dist < 1.5
        powerLines_pro(i) = [];
    else
        i = i + 1;
    end
end

% Calculate the length of the power line
i = 1;
sumdist = 0;
while i <= size(powerLines_pro,2)
    dist = getDist(powerLines_pro(i).Location);
    i = i + 1;
    sumdist = sumdist + dist;
end
sumdist

% Insert points to sparse power lines
for i = 1:size(powerLines_pro,2)
    powerLines_pro_new(i).Location = insert_3D(powerLines_pro(i).Location, 0.1);
end
% save('C:\Users\Lily\Desktop\PLM','powerLines_pro_new')

% powerLines_pro_new = powerLines_pro;
% visualization
powerLines_new = [];
colors = [];
for i = 1:size(powerLines_pro_new,2)
    powerLines_new = [powerLines_new; powerLines_pro_new(i).Location];
    temp = repmat(rand(1,3),size(powerLines_pro_new(i).Location,1),1);
    colors = [colors; temp];
    PLMs(i).Location = powerLines_pro_new(i).Location;
    PLMs(i).Color = temp;
end
% ptplm = pointCloud(powerLines_new,'Color',colors)
figure
% pcshow(ptplm.Location(1:1:end,:),ptplm.Color(1:1:end,:))
pcshow(powerLines_new,colors)
title('Power line modeling')
view(60,25)
print(gcf,'-dpng','-r300', 'f7_Power line model.png')













































