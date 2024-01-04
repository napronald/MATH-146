clear all; clc;

T = [ 0   0   1/2
      1/2 1/3 0
      1/2 2/3 1/2]

lamb1 = 1/2;
lamb2 = 1/3;
lamb3 = 1/4;

A = zeros(4,3);
A(1,1) = -lamb1;
A(2,2) = -lamb2;
A(3,3) = -lamb3;

A(1,2) = T(4)*lamb2;
A(1,3) = T(7)*lamb3;

A(2,1) = T(2)*lamb1;
A(2,3) = T(8)*lamb3;

A(3,1) = T(3)*lamb1;
A(3,2) = T(6)*lamb2;

b = [0;0;0;1]
A(4,:) = 1;
A

A'

A'*A \ A'*b

