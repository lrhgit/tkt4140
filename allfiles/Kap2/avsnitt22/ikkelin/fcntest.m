function yd = fcntest(x,y)
%--- Brukes av programmettest
% Setter a = 6
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = 6*y(1)^2;