%% show 
function result = showbreaks(A)
A = A(1:10:size(A,1),:);
A_shift = A-mean(A);
[eValue,eVector,angle] = eigenDV(A_shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
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