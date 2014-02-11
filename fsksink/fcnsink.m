function yd = fcnsink(x,y)
% Called from  SINK and FCNPHI 
yd = zeros(3,1);
yd(1) = y(2);
yd(2) = y(3);
yd(3) = y(2)^2 - 1;