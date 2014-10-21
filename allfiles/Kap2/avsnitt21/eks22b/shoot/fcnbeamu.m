function dy = fcnbeamu(x,y)
% Used by program udis
dy = zeros(size(y));
dy(1) = y(2);
dy(2) = -(1 + x^2)*(y(1) + (1.0 - x^2)/2) ; 
