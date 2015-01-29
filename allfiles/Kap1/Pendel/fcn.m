function dydt = fcn(t,y)
% Kalles av programmet Pendel
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -sin(y(1));

