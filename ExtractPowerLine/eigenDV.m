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
cross_vector = cross(N,eVector); % The Z component is regular in the positive Z direction
if(cross_vector(3)<0)
    eVector = -eVector;
end
angle = acos(sum(N.*eVector)/sqrt(sum(N.^2))*sqrt(sum(eVector.^2)))*180.0/pi;
end