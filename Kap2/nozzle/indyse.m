% ======================= program indyse ==========================
% Løser dyseproblemet som et intialproblem
% Ligning : f'''(x)+ f*f''(x) + (f'(x))^2 = 0
%            f(0) = 0, f'(0) = 2 ,f''(0)= 0
% Bruker funksjonen fcndyse
% =============================================================
clear
xinf = 15.0;
xspan = [0.0 xinf];
f0 = [0 ; 0.1; 0];
options = odeset('RelTol',1.0e-5);
[x,f] = ode45(@fcndyse,xspan,f0,options);
plot(x,f(:,1),x,f(:,2),x,f(:,3));
grid
shg
%---- Compute a table of w and w'
% zspan = [0:0.2:zinf];
% w0 = [1 ; s];
% options = odeset('RelTol',1.0e-6);
% [z,w] = ode45(@fcn51,zspan,w0,options,alpha);
% fprintf('\n           z         w           w''\n\n');
% fprintf(' %12.2f %12.5f %12.5f\n',[z w]');
