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
