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