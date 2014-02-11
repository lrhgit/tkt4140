function yd = fcnjh(x,y)
% Used by jhamel and plotjh
% solving Jeffrey-Hamel's equation
global Re alpha
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -2*alpha*((Re*y(1) + 2*alpha)*y(2));