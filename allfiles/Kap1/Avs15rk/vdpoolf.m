function dydt = vdpoolf(t,y,mu )
%UNTITLED2 Summary of this function goes here
% 2x2 system for a van der Pool oscillator.

dydt = zeros(length(y),1);
dydt(1) = y(2);
dydt(2) = mu*(1-y(1)^2)*y(2)-y(1);
end

