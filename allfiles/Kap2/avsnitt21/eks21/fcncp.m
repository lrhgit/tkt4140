function dy = fcncp(x,y)
% Used by cpskyt
dy = zeros(size(y));
global P ;
dy(1) = y(2); % 1. equation
dy(2) = - P;  % 2. equation

