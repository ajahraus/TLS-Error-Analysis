function [Rotation, trans] = HornsRegistration(child, parent)
% [Rotation, trans] = HornsRegistration(child, parent)
% This function returns the rotation matrix R and translation vector t,
% such that b = a*R + t; where b is the Parent coordinates [x,y,z], a is 
% the child coordinates, and each row of a corresponds to the same point in
% b.
% a = [x1, y1, z1; x2, y2, z2; ... x_n, y_n, z_n], nx3 matrix
% b has the same format as a;
% R = R3(kappa)R2(phi)R1(omega), a 3x3 matrix;
% t = [x_t, y_t, z_t];
% IMPORTANT: If a and b are column vectors rather than row vectors, then
% the expression becomes R'*a + t' = b

%
centroid1 = sum(child)/size(child,1);
centroid2 = sum(parent)/size(parent,1);

a(:,1) = child(:,1) - centroid1(1);
a(:,2) = child(:,2) - centroid1(2);
a(:,3) = child(:,3) - centroid1(3);

b(:,1) = parent(:,1) - centroid2(1);
b(:,2) = parent(:,2) - centroid2(2);
b(:,3) = parent(:,3) - centroid2(3);

%
sxx = sum(a(:,1).*b(:,1));
sxy = sum(a(:,1).*b(:,2));
sxz = sum(a(:,1).*b(:,3));

syx = sum(a(:,2).*b(:,1));
syy = sum(a(:,2).*b(:,2));
syz = sum(a(:,2).*b(:,3));

szx = sum(a(:,3).*b(:,1));
szy = sum(a(:,3).*b(:,2));
szz = sum(a(:,3).*b(:,3));

N = zeros(4);
N(1,1) = sxx+syy+szz; N(2,2) = sxx - syy - szz;
N(3,3) = syy - sxx - szz; N(4,4) = szz - sxx - syy;
N(2,1) = syz - szy; N(1,2) = N(2,1); N(3,1) = szx - sxz; N(1,3) = N(3,1);
N(4,1) = sxy - syx; N(1,4) = N(4,1); N(2,3) = sxy + syx; N(3,2) = N(2,3);
N(2,4) = szx + sxz; N(4,2) = N(2,4); N(3,4) = syz + szy; N(4,3) = N(3,4);

[~ , idx ]= max(eig(N));
[V,~] = eig(N);
v = V(:,idx);

M = zeros(3);
M(1,:) = [v(1)^2 + v(2)^2 - v(3)^2 - v(4)^2, 2*(v(2)*v(3) - v(1)*v(4)), 2*(v(2)*v(4) + v(1)*v(3))];
M(2,:) = [2*(v(2)*v(3) + v(1)*v(4)),v(1)^2 - v(2)^2 + v(3)^2 - v(4)^2, 2*(v(3)*v(4) - v(1)*v(2))];
M(3,:) = [2*(v(2)*v(4) - v(1)*v(3)),2*(v(3)*v(4) + v(1)*v(2)),v(1)^2 - v(2)^2 - v(3)^2 + v(4)^2];

Rotation = M';
trans = centroid2-(centroid1*Rotation);