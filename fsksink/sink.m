% === Program SINK
% Solution of the Falkner-Skan equation using 
% a shooting technique for the case with beta -> infinity
% This is the only case having an analytical solution.
%  The equation is given by :
%  g'''(ksi) + [1 - (g'(ksi))^2] = 0
%  g(0) = 2*(sqrt(3) - 1), g'(0) = -4/3, g''(0) = 14/9 and g'(ksinf) = 1
% In the program vi put g = y(1), g' =  y(2) and g'' = y(3). 
% and we use x instead of ksi
% Using FZERO instead of the secant-method.
% Boundary conditions used :
% g(0) = 2*(sqrt(3) - 1) , g' = -4/3, g'(ksinf) = 1
% s = g''(0)
clear
global x1;
x1 = 6.0; %ksi-inf
s0 = 1.5; % guessing
tols = 1.0e-8;
xspan = [0 x1];
%options = optimset('TolX',tols);
options = optimset('Display','iter','TolX',tols);
s = fzero(@fcnphi,s0,options);
fprintf('Computed value of g"(0) = %12.5e \n',s);

% Table of g, g' og g'' 
xspan = [0:0.25:x1];
g0 = 2*(sqrt(3) - 1);
y0 = [g0 ; -4/3 ; s];
tol = 1.0e-9;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol*1.0e-4]);
[x,y] = ode15s(@fcnsink,xspan,y0,options);
fprintf('\n         ksi        g          g''        g"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[x y]');
