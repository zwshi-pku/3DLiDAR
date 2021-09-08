%% Fit wires, and insert points at sparse parts with fixed resolution to form new wire data points.
function xyzs_new = insert_3D(cluster_raw, resolution)
cluster_shift = cluster_raw-mean(cluster_raw);
if size(cluster_shift,1)<3
    xyzs_new = cluster_raw;
    return;
end

[eValue,eVector,angle] = eigenDV(cluster_shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
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