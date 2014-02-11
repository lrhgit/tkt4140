function dydt = fcndpend(t,y)
% Brukes av programmet dpendel
global g m1 m2 l1 l2;
dydt = zeros(size(y));
dydt(1) = y(2);
dydt(2) = m2*l2*sin(y(3) - y(1))*y(4)^2 - g*(m1+m2)*sin(y(1));
dydt(3) = y(4);
dydt(4) = -l1*sin(y(3) - y(1))*y(2)^2 - g*sin(y(3));