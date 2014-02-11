function yd = fcnbla(x,y)
% Used by blaplot and blasec
yd = zeros(3,1);
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -y(1)*y(3);