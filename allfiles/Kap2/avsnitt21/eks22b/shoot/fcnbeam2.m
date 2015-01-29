function dy = fcnbeam2(x,y)
% Used by program momentdis
global n;
dy = zeros(size(y));
dy(1) = y(2);
dy(2) = -(y(1)*(1 + x^n) + 1.0) ; 
