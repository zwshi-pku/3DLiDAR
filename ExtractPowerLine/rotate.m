%% Rotate a certain number of degrees counterclockwise about the Z axis
function [cloudR R]= rotate(cloud, alpha)
R = [cos(alpha) -sin(alpha) 0;
    sin(alpha) cos(alpha) 0;
    0 0 1];
cloudR = cloud*R';
end