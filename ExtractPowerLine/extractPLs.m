function isPLIndex = extractPLs(nonGroundPoints,radius,angleThr,LThr)
[normals, Ls] = getPCA(nonGroundPoints,radius);
angle = acosd(normals(:,3)./sqrt(sum(normals.^2,2)));
isPLIndex = abs(angle - 90) < angleThr & Ls>LThr;

% normals = pcnormals(pointCloud(nonGroundPoints),30);
% figure
% pcshow(nonGroundPoints)
% title('Adjusted Normals of Point Cloud')
% hold on
% quiver3(nonGroundPoints(:,1),nonGroundPoints(:,2),nonGroundPoints(:,3), normals(:,1), normals(:,1), normals(:,1));
% hold off
end