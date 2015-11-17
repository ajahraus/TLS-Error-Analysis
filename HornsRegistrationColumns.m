function [Rotation, trans] = HornsRegistrationColumns(child, parent)
% [Rotation, trans] = HornsRegistration(child, parent)
% This function returns the rotation matrix R and translation vector t,
% such that b = R*a + t; where b is the Parent coordinates [x,y,z], a is 
% the child coordinates, and each column of a corresponds to the same point in
% b.
% a = [x1, x2, x3...;y1, y2, y3; z1, z2, z3;], 3xn matrix
% b has the same format as a;
% R = R3(kappa)R2(phi)R1(omega), a 3x3 matrix;
% t = [x_t; y_t; z_t];

%
centroid1 = mean(child,2);
centroid2 = mean(parent,2);

a = child-repmat(centroid1,1,size(child,2));
b = parent-repmat(centroid2,1,size(parent,2));

%
sxx = a(1,:)*b(1,:)';
sxy = a(1,:)*b(2,:)';
sxz = a(1,:)*b(3,:)';

syx = a(2,:)*b(1,:)';
syy = a(2,:)*b(2,:)';
syz = a(2,:)*b(3,:)';

szx = a(3,:)*b(1,:)';
szy = a(3,:)*b(2,:)';
szz = a(3,:)*b(3,:)';

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

Rotation = M;
trans = mean(parent - Rotation*child,2);