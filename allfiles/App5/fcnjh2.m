function yd = fcnjh2(x,y)
% Used by jhfigur when solving
% Jeffrey-Hamel's equation for the case 
% when alpha^2 << Re*alpha
global Rea 
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = - 2*Rea*y(1)*y(2);