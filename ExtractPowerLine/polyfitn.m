%% The function order is obtained by automatic fitting：[p,s,m,n] = polyfitn(x,y,10,0.01)
function [p,n] = polyfitn(x,y,top_n,thr_precise)
% for i=1:top_n
%     p = polyfit(x,y,i);
%     f = polyval(p,x);     %Calculate the value of the fitting function at X.
%     if sum((f-y).^2)<thr_precise
%         n = i;
%         break;
%     end
% end
% if i==top_n
%     disp('The fitting order is not obtained');
%     n = 2;
% end
n = 2;
p = polyfit(x,y,n);
end
