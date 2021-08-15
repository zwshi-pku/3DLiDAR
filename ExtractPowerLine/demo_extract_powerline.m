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
ptpl = pointCloud(powerLines_new,'Color',colors)
figure
pcshow(ptpl)
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
ptplm = pointCloud(powerLines_new,'Color',colors)
figure
pcshow(ptplm.Location(1:1:end,:),ptplm.Color(1:1:end,:))
title('Power line modeling')
view(60,25)
print(gcf,'-dpng','-r300', 'f7_Power line model.png')















%% Calculate power line length
function dist = getDist(powerLines)
if size(powerLines,1)<3
    dist = 0;
    return;
end
shift = powerLines-mean(powerLines);
[eValue,eVector,angle] = eigenDV(shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
rotated = rotate(shift, -angle*pi/180.0);
shift_x = rotated(:,1);
dist = max(shift_x)-min(shift_x);
end

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

% 删除被标记的数据
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







%% Rotate a certain number of degrees counterclockwise about the Z axis
function [cloudR R]= rotate(cloud, alpha)
R = [cos(alpha) -sin(alpha) 0;
    sin(alpha) cos(alpha) 0;
    0 0 1];
cloudR = cloud*R';
end

%% The function order is obtained by automatic fitting：[p,s,m,n] = polyfitn(x,y,10,0.01)
function [p,n] = polyfitn(x,y,top_n,thr_precise)
% for i=1:top_n
%     p = polyfit(x,y,i);
%     f = polyval(p,x);     %计算拟合函数在x处的值。
%     if sum((f-y).^2)<thr_precise
%         n = i;
%         break;
%     end
% end
% if i==top_n
%     disp('未得到拟合阶数');
%     n = 2;
% end
n = 2;
p = polyfit(x,y,n);
end

%% Get the maximum eigenvalue and its corresponding eigenvector, and return the included angle with the X axis
function [eValue,eVector,angle] = eigenDV(xyzs)
covm = cov(xyzs);
[V,D] = eig(covm);
[D_sort,index] = sort(diag(D),'descend');
V_sort = V(:,index);

eVector = V_sort(:,1)';
eValue = D_sort(1);
X = [1 0 0];
Y = [0 1 0];
Z = [0 0 1];
N = X;
cross_vector = cross(N,eVector); % z分量正则朝向z正方向
if(cross_vector(3)<0)
    eVector = -eVector;
end
angle = acos(sum(N.*eVector)/sqrt(sum(N.^2))*sqrt(sum(eVector.^2)))*180.0/pi;
end

%% Add data point
function [points, insert_pnts]= insert(xys, p, resolution)
insert_pnts = [];
minValue = min(xys(:,1));
maxValue = max(xys(:,1));
xs_gap = linspace(maxValue, minValue, floor(abs(maxValue-minValue)/resolution));
if isa(p,'cfit')
    ys_gap = p(xs_gap);
    points = [xs_gap' ys_gap];
else
    ys_gap = polyval(p,xs_gap);
    points = [xs_gap' ys_gap'];
end

insert_pnts = points;
return;

% insert_pnts = [];
% i = 1;
% num = size(xys,1)-1;
% while i<=num
%     dist = sqrt(sum((xys(i,:)-xys(i+1,:)).^2));
%     n = floor((dist/resolution));
%     if n>0
%         xs = linspace(xys(i,1),xys(i+1,1),n + 2);
%         xs(end) = [];
%         xs(1) = [];
%         ys = polyval(p,xs);
%         xys = [xys(1:i,:);[xs' ys'];xys(i+1:end,:)];
%         insert_pnts = [insert_pnts; [xs' ys']];
%     end
%     i = i + 1;
% end
% points = xys;
end

%% Add data point
function [points, insert_pnts]= insert3(xyzs, p, resolution)
insert_pnts = [];
i = 1;
num = size(xyzs,1)-1;
while i<=num
    dist = sqrt(sum((xyzs(i,:)-xyzs(i+1,:)).^2));
    n = floor((dist/resolution));
    if n>0
        xs = linspace(xyzs(i,1),xyzs(i+1,1),n + 2);
        xs(end) = [];
        xs(1) = [];
        zs = polyval(p,xs);
        inserted = [xs' zeros(size(xs,2),1) zs'];
        xyzs = [xyzs(1:i,:);inserted;xyzs(i+1:end,:)];
        insert_pnts = [insert_pnts; inserted];
    end
    i = i + 1;
end
points = xyzs;
end

%% Fit wires, and insert points at sparse parts with fixed resolution to form new wire data points.
function xyzs_new = insert_3D(cluster_raw, resolution)
cluster_shift = cluster_raw-mean(cluster_raw);
if size(cluster_shift,1)<3
    xyzs_new = cluster_raw;
    return;
end

[eValue,eVector,angle] = eigenDV(cluster_shift);

% 绕z轴顺时针旋转一定度数，数据B在X-Z直角坐标内
rotated = rotate(cluster_shift, -angle*pi/180.0);
x = rotated(:,1);
y = rotated(:,2);
z = rotated(:,3);
p = polyfit(x,z,2);
% p = catenary(x,z);
[xzs_new xzs_inserted] = insert([x z], p, resolution);
if isempty(xzs_inserted)
    xyzs_new = cluster_raw;
    return;
end
% xyzs_new = [x y z; xzs_inserted(:,1) zeros(size(xzs_inserted,1),1) xzs_inserted(:,2)];
xyzs_new = [xzs_inserted(:,1) zeros(size(xzs_inserted,1),1) xzs_inserted(:,2)];

% [x_new,ind] = sort(xyzs_new(:,1),'descend');
% xyzs_new = xyzs_new(ind,:);
xyzs_new = rotate(xyzs_new, angle*pi/180.0)+mean(cluster_raw);
end

%% Merge over split powerline points
function ids = findMerge(powerLines, index, begin)
ids = [];
A = powerLines(index(begin)).Location;
if size(A,1)<10
   return; 
end
meanAz = mean(A(:,3));
% 
A_shift = A-mean(A);
[eValue,eVector,angle] = eigenDV(A_shift);

% 绕z轴顺时针旋转一定度数，数据B在X-Z直角坐标内
ARotated = rotate(A_shift, -angle*pi/180.0);
A_shift_x = ARotated(:,1);
A_shift_z = ARotated(:,3);
p = polyfit(A_shift_x,A_shift_z,2);
% p = catenary(A_shift_x,A_shift_z);
pl = polyfit(A(:,1),A(:,2),1);

for i = begin+1:size(index,1)
    if powerLines(index(i)).Label == 1
        continue;
    end
    
    B = powerLines(index(i)).Location;
    meanBz = mean(B(:,3));
    PLB_shift = B-mean(A);
    BRotated = rotate(PLB_shift, -angle*pi/180.0);
    B_shift_x = BRotated(:,1);
    B_shift_y = BRotated(:,2);
    B_shift_z = BRotated(:,3);
%     meanBy = abs(mean(B_shift_y))
    
    if isa(p,'cfit')
        B_Pz = p(B_shift_x);
    else
        B_Pz = polyval(p,B_shift_x);
    end
    B_Py = polyval(pl,B(:,1));
    deltaBy = abs(mean(B_Py - B(:,2)));

%     meanBd = abs(mean(B_Pz - B_shift_z));
    deltaBz = abs(meanAz-meanBz);
    meanBd = max(abs(B_Pz - B_shift_z));
    if deltaBy < 0.5 & meanBd < 0.2 & deltaBz < 4
        result = 1;
    else
        result = 0;
    end
    if 0
        figure(3)
        pcshow(ARotated)
        hold on
        pcshow(BRotated)
        close 3;        
    end
    clc

    if result == 1
        ids = [ids; index(i)];
        A_shift_x = [A_shift_x; B_shift_x];
        A_shift_z = [A_shift_z;B_shift_z];
        p = polyfit(A_shift_x,A_shift_z,2);
    end
end
end

function result = isMerge(PLA, PLB)
PLA_shift = PLA-mean(PLA);
if size(PLA_shift,1)<3
    result = 0;
    return;
end

[eValue,eVector,angle] = eigenDV(PLA_shift);

% 绕z轴顺时针旋转一定度数，数据B在X-Z直角坐标内
rotated = rotate(PLA_shift, -angle*pi/180.0);
PLA_shift_x = rotated(:,1);
PLA_shift_y = rotated(:,2);
PLA_shift_z = rotated(:,3);
meanPLAy = mean(PLA_shift_y);

PLB_shift = PLB-mean(PLA);
rotated = rotate(PLB_shift, -angle*pi/180.0);
PLB_shift_x = rotated(:,1);
PLB_shift_y = rotated(:,2);
PLB_shift_z = rotated(:,3);
meanPLBy = mean(PLB_shift_y);

% p = polyfit(PLA_shift_x,PLA_shift_z,2);
% PLB_Pz = polyval(p,PLB_shift_x);
p = catenary(PLA_shift_x,PLA_shift_z);
PLB_Pz = p(PLB_shift_x);
mean_dist = abs(mean(PLB_Pz - PLB_shift_z));

if abs(meanPLBy) < 0.1 & mean_dist < 0.1
    result = 1;
else
    result = 0;
end
end

%% Catenary: y=a+ccosh((x-b)/c)
function cfun = catenary(xs,ys)
syms x;
f = fittype('a+c*cosh((x-b)/c)','independent','x','coefficients',{'a','c','b'});
cfun = fit(xs,ys,f,'StartPoint', [-780.6,780.4,-1.184])
end

%% show 
function result = showbreaks(A)
A = A(1:10:size(A,1),:);
A_shift = A-mean(A);
[eValue,eVector,angle] = eigenDV(A_shift);

% 绕z轴顺时针旋转一定度数，数据B在X-Z直角坐标内
rotated = rotate(A_shift, -angle*pi/180.0);
[temp, ind] = sort(rotated(:,1),'ascend');
rotated = rotated(ind,:);

C1 = rotated(1:floor(size(rotated,1)*0.5),:);
C2 = rotated(floor(size(rotated,1)*2/3):end,:);
C3 = rotated(floor(size(rotated,1)*0.5):floor(size(rotated,1)*2/3),:);
figure
set (gcf,'position',[200,200,700,700] );
hold on
box on
plot3(C1(:,1),C1(:,2),C1(:,3),'bo');
plot3(C2(:,1),C2(:,2),C2(:,3),'co');
% plot3(C3(:,1),C3(:,2),C3(:,3),'go');
pbaspect([20 1 2])
xlabel('\itX^{''}','FontSize',10,'FontName','times new Roman')
ylabel('\itY^{''}','FontSize',10,'FontName','times new Roman')
zlabel('\itZ^{''}','FontSize',10,'FontName','times new Roman')
print(gcf,'-dpdf','-r300','C:\Users\Lily\Desktop\linebreaks.pdf');
cf1 = catenary(C1(:,1),C1(:,3));
cf_x = rotated(:,1);
cf_z = cf1(cf_x);
cf_y = zeros(size(cf_x,1),1);
plot3(cf_x,cf_y,cf_z,'b-','linewidth',1);

C = [C1;C2];
cf2 = catenary(C(:,1),C(:,3));
c_x = rotated(:,1);
c_z = cf2(c_x);
c_y = zeros(size(c_x,1),1);
hold on
plot3(c_x,c_y,c_z,'r-','linewidth',1);
print(gcf,'-dpdf','-r300','C:\Users\Lily\Desktop\linefit.pdf');

minx = min(C3(:,1));
maxx = max(C3(:,1));
gap_x = linspace(minx,maxx,floor(abs(maxx-minx)/1.0));
gap_y = zeros(size(gap_x,2),1);
gap_z = cf2(gap_x');
hold on
plot3(gap_x,gap_y,gap_z,'ro');
print(gcf,'-dpdf','-r300','C:\Users\Lily\Desktop\points_gap.pdf');
end

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

% 判断聚类是否可以拟合悬链模型
% result = 0模型拒绝，result = 1模型采纳。
function mse = isFit(A)
mse = 10;
if size(A,1)<10
    return;
end

A_shift = A-mean(A);
[eValue,eVector,angle] = eigenDV(A_shift);

% 绕z轴顺时针旋转一定度数，数据B在X-Z直角坐标内
ARotated = rotate(A_shift, -angle*pi/180.0);
A_shift_x = ARotated(:,1);
A_shift_z = ARotated(:,3);
p = polyfit(A_shift_x,A_shift_z,2);
% p = catenary(A_shift_x,A_shift_z);
A_new_z = polyval(p,A_shift_x);
mse = rmse(A_shift_z,A_new_z);
end
