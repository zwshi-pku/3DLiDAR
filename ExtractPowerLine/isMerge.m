function result = isMerge(PLA, PLB)
PLA_shift = PLA-mean(PLA);
if size(PLA_shift,1)<3
    result = 0;
    return;
end

[eValue,eVector,angle] = eigenDV(PLA_shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
rotated = rotate(PLA_shift, -angle*pi/180.0);
PLA_shift_x = rotated(:,1);
PLA_shift_y = rotated(:,2);
PLA_shift_z = rotated(:,3);
meanPLAy = mean(PLA_shift_y);

PLB_shift = PLB-mean(PLA);
rotated = rotate(PLB_shift, -angle*pi/180.0);
PLB_shift_x = rotated(:,1);
PLB_shift_y = rotated(:,2);
PLB_shift_z = rotated(:,3);
meanPLBy = mean(PLB_shift_y);

% p = polyfit(PLA_shift_x,PLA_shift_z,2);
% PLB_Pz = polyval(p,PLB_shift_x);
p = catenary(PLA_shift_x,PLA_shift_z);
PLB_Pz = p(PLB_shift_x);
mean_dist = abs(mean(PLB_Pz - PLB_shift_z));

if abs(meanPLBy) < 0.1 & mean_dist < 0.1
    result = 1;
else
    result = 0;
end
end