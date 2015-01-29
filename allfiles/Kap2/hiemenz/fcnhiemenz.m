function yd = fcnhiemenz(x,y)
% Brukes av hiemenzplot
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -(y(1)*y(3) + (1 - y(2)^2));