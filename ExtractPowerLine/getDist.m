%% Calculate power line length
function dist = getDist(powerLines)
if size(powerLines,1)<3
    dist = 0;
    return;
end
shift = powerLines-mean(powerLines);
[eValue,eVector,angle] = eigenDV(shift);

% Rotate clockwise around the Z axis for a certain number of degrees, and the data B is in the x-z rectangular coordinates
rotated = rotate(shift, -angle*pi/180.0);
shift_x = rotated(:,1);
dist = max(shift_x)-min(shift_x);
end