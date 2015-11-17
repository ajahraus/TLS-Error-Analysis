function R = rotzd(a)

R = eye(3);

R(1:2,1:2) = [cosd(a), -sind(a); sind(a), cosd(a)];
