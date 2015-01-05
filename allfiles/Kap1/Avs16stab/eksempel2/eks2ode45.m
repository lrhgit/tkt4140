function eks2ode45
% Løser eksempel 2 i avsnitt 1.6.2 i kompendiet
% Eksempel på et stivt system
% Bruker ode45 istedenfor ode15s som i eks2ode15s
% 
clear
options = odeset('RelTol',1.0e-5,'AbsTol', 1.0e-8,'Refine',1);
                 
tspan = [0 0.12516];
   y0 = [0.0 0.0 ];
   [t,y] = ode45(@f,tspan,y0,options);
   fprintf('        t           y              v \n\n');
   fprintf('  %12.4e  %12.5e  %12.5e\n',[t y]');
 length(t)
function dydt = f(t,y)   
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -129600.0*y(2) + 98696.0*y(1)+ 9869.6;
% -----------------------------------------------

