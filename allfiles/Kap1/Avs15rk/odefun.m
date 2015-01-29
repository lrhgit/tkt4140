function dydt = odefun(t,y,g,a)
% Used by rkneks13
% Falling sphere with constant Cd
dydt = zeros(length(y),1);
dydt(1) = y(2);
dydt(2) = g - a*y(2)^2;
