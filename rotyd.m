function R = rotyd(a)

R = eye(3);

R(1,1) = cosd(a);
R(1,3) = sind(a);
R(3,3) = cosd(a);
R(3,1) = -sind(a);
