%===================== eulv =============================
% The program computes the velocity of a golf ball
% falling vertically in air. The equation 
% of motion is integrated using the forward Euler- method
%========================================================

clear
% Using Cd = 0.4;
g = 9.81    ;  % Gravity [N/kg]
alfa = 7.0e-3;
tend = 10.0  ; % Max. time [s]
dt = 0.5; % Timestep [s]
steps = round(tend/dt) + 1;
v = zeros(steps,1); t = v;  % allocate space
v(1)= 0.0 ; t(1) = 0.0;
for n = 1:steps - 1
   t(n+1) = n*dt;
   v(n+1) = v(n) + dt*(g - alfa*v(n)^2);
end
% Analytical solution 
k1 = sqrt(g/alfa); k2 = sqrt(alfa*g);
va = k1*tanh(k2*t);
fprintf('       t(s)       v(m/s)     va(m/s) \n\n');
fprintf(' %10.2f  %10.3f %10.3f \n',[t v va]');
clf
plot(t,v,'k',t,va,'k-.');
grid on
FS = 'FontSize';
xlabel('t(s)',FS,14)
ylabel('v, v_{a}(m/s)',FS,12)
legend('v','v_{a}')
title('Euler''s method with \Deltat = 0.5',FS,14)
shg
