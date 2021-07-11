function [K] = SpringStiffnessMatrix(V,k)

%V = 2x3 matrix of the global coords of the four nodes
%k = 6x1 vector of the six spring constants

L = sqrt(sum((V(1,:)-V(2,:)).^2));

Cx = (V(1,1) - V(2,1))/L;
Cy = (V(1,2) - V(2,2))/L;
Cz = (V(1,3) - V(2,3))/L;

A = [Cx^2 Cx*Cy Cx*Cz; Cy*Cx Cy^2 Cy*Cz; Cz*Cx Cz*Cy Cz^2;];

K = k*[A -A; -A A;];

end