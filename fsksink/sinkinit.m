% === Program SINKINIT
% Numerical test of the Falkner-Skan equation 
% for the case with beta -> infinity
% This is the only case having an analytical solution.
%  The equation is given by :
%  g'''(ksi) + [1 - (g'(ksi))^2] = 0
%  g(0) = 2*(sqrt(3) - 1), g'(0) = -4/3, g''(0) = 14/9 and g'(ksinf) = 1
% In the program vi put g = y(1), g' =  y(2) and g'' = y(3). 
% and vi use x instead of ksi
%
clear
x0 = 0; %start
x1 = 7.5; %ksinf
xspan = [x0:0.25:x1];
% Startverdier
y01 = 2.0*(sqrt(3) - 1.0); y02 = -4/3; y03 = 14/9;
y0 = [y01 ; y02 ; y03];
tol = 1.0e-10;
options = odeset('RelTol',tol,'AbsTol',[tol tol tol*1.0e-3] );
[x,y] = ode45('fcnsink',xspan,y0,options);

fprintf('\n         ksi        g          g''        g"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[x y]');

% Analytiske verdier
v1 = x; % initialiser
v1 = tanh((x/sqrt(2)) + atanh(sqrt(2)/3));
g = x + 2*sqrt(3) - 3*sqrt(2)*v1;
dg = -2 + 3*v1.*v1;
d2g = sqrt(2)*3*v1.*(1-v1.*v1);
fprintf('\n Analytiske verdier \n')
fprintf('\n         ksi        g          g''        g"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[x g dg d2g]');
% Plotting
clf
plot(x,g,x,dg,x,d2g)
grid on
xlabel('\xi','FontSize',14,'FontWeight','Bold')
ylabel('g, g'' , g"','Fontsize',14)
title('Fsk-Sink','Fontsize',14)
%legend('f''','f"')
