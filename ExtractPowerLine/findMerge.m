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

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
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