function eks1ode15s
% Løser eksempel 1 i avsnitt 1.6.1 i kompendiet
% Eksempel på et stivt system.
% Bruker ode15s som løser
% Må kjøres under versjon 7.x av Matlab
clear
c = 100;
l1 = -2*c; l = -1;
options = odeset('RelTol',1.0e-5,'AbsTol', 1.0e-8,'Jacobian',@Jac);
                 
tspan = [0 20];
   x0 = [1.0 0.0 ];
   [t,x] = ode15s(@f,tspan,x0,options);
   fprintf('        t           x              v \n\n');
   fprintf('  %12.4e  %12.5e  %12.5e\n',[t x]');
   % length(t)
 
function dxdt = f(t,x)   
dxdt = zeros(2,1);
dxdt(1) = x(2);
dxdt(2) = l1*x(2) + l*x(1);
end
% ------------------------------------------------
function dfdx = Jac(t,x)
dfdx = [0 1 ; l l1];
end
end

