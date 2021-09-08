function [normals, Ls] = getPCA(xyz,r)
idxNN = rangesearch(xyz, xyz, r);
normals   = nan(size(xyz,1),3);
Ls = nan(size(xyz,1),1);

for i=1:size(xyz,1)
    if numel(idxNN{i})<3
        continue;
    end
    XNN = xyz(idxNN{i},:);
    
    covm = cov(XNN);
    [V, lambda] = pcacov(covm);
    normals(i,:) = V(:,1)';
    lmbda1 = lambda(1);
    lmbda2 = lambda(2);
    lmbda3 = lambda(3);
    
%     [V,D] = eig(covm);
%     [D_S,index] = sort(diag(D),'descend');
%     V_S = V(:,index);
%     
%     lmbda1 = D_S(1);
%     lmbda2 = D_S(2);
%     lmbda3 = D_S(3);
    
    L = (lmbda1-lmbda2)/lmbda1;
    P = (lmbda2-lmbda3)/lmbda1;
    S = lmbda3/lmbda1;
    Ls(i,1) = L;
end
end