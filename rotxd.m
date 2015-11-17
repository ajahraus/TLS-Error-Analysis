function R = rotxd(a)

R = eye(3);

R(2:3,2:3) = [cosd(a), -sind(a); sind(a), cosd(a)];
