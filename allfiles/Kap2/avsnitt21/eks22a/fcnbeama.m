function dy = fcnbeama(x,y)
% Used by program beamcola
dy = zeros(size(y));
global theta2 velg ;
dy(1) = y(2);
dy(2) = - theta2*y(1); % 2. equation
if( velg ==1)
   dy(2) = dy(2) - theta2*(1.0 - x^2)/2; % 1. equation
end
