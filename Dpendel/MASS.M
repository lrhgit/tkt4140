function M = mass(t,y)
% Massematrisa for program dpendel
global g m1 m2 l1 l2;
n = length(y);
M = zeros(n,n);
M(1,1) = 1.0;
M(2,2) = l1*(m1 + m2);
M(2,4) = m2*l2*cos(y(3) - y(1));
M(3,3) = 1.0;
M(4,2) = l1*cos(y(3) - y(1));
M(4,4) = l2;