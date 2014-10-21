function dy = fcnknekk(x,y)
% Used by program beamcol
dy = zeros(size(y));
global theta2 velg ;
dy(1) = y(2);
dy(2) = - theta2*4*y(1); % 2. equation
if( velg ==1)
   dy(2) = dy(2) - theta2*2*x*(1.0 - x); % 1. equation
end
