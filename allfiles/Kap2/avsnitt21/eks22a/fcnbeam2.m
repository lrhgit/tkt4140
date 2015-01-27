function dy = fcnbeam2(x,y)
% Used by program beamcolv2
dy = zeros(size(y));
global theta2;
dy(1) = y(2);
dy(2) = -theta2*(y(1) + (1.0 - x^2)*0.5); 
