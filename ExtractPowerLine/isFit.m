% Judge whether clustering can fit the catenary model
% result = 0 reject，result = 1 accept
function mse = isFit(A)
mse = 10;
if size(A,1)<10
    return;
end

A_shift = A-mean(A);
[eValue,eVector,angle] = eigenDV(A_shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
ARotated = rotate(A_shift, -angle*pi/180.0);
A_shift_x = ARotated(:,1);
A_shift_z = ARotated(:,3);
p = polyfit(A_shift_x,A_shift_z,2);
% p = catenary(A_shift_x,A_shift_z);
A_new_z = polyval(p,A_shift_x);
mse = rmse(A_shift_z,A_new_z);
end