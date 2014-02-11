function eks2ode15s
% Løser eksempel 2 i avsnitt 1.6.2 i kompendiet
% Eksempel på et stivt system med hendelse
% Bruker ode15s som løser
clear
options = odeset('RelTol',1.0e-5,'AbsTol', 1.0e-8,'Jacobian',@Jac,...
                 'Events',@events);
tspan = [0 0.15];
   y0 = [0.0 0.0 ];
   [t,y,te,ye,ie] = ode15s(@f,tspan,y0,options);
   fprintf('        t           y              v \n\n');
   fprintf('  %12.4e  %12.5e  %12.5e\n',[t y]');
   % length(t)
 
function dydt = f(t,y)   
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -129600.0*y(2) + 98696.0*y(1)+ 9869.6;
% -----------------------------------------------
function [value,isterminal,direction] = events(t,y)
% Finn tiden når partikkelen treffer ytre sylinder
value = [0.01 - y(1); y(2)];
isterminal = [1; 0];
direction = [0 ; 0];
% ------------------------------------------------
function dfdy = Jac(t,y)
dfdy = [0 1 ; 98696.0 -129600.0];
% -----------------------------------------------

