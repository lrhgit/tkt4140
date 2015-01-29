function a = accx(x,y )
%Gravity force in direction x
%   Detailed explanation goes here
%
global MG
r=sqrt(x^2 + y^2);
a = -MG*x/r^3;
end

