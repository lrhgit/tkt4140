%  ============== Program blasec ===============
%  Solution of Blasius equation using a shooting
%  technique.
% (See  example 2.6 in the compendium, ver. 2010))
%
%  Blasius equation is given by:
%  f'''(eta) + f(eta)*f''(eta) = 0
%  f(0) = 0, f'(0) = 0, f'(etainf)=1
% In the program we put y(1)= f, y(2) = f'
% and y(3) = f'' and using  x instead of  eta
% We put etainf = 5.75 (x1 in the program)
% s = f''(0) = y(3)(0)
%
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20); %Default values for plotting.

x0 = 0; %start
x1 = 5.75; %etainf
xspan = [x0 x1];


%% From the program blaplot we have found
%  the starting values s0 and s1 for f''(0)

s0 = 0.3;
s1 = 0.5;
% Compute fi0
y0 = [0.0 0.0  s0];
[x,y] = ode45(@blasius,xspan,y0);
fi0 = y(end,2) - 1;

%% Initial values for the iteration
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
% Headline of the table
fprintf('        itr.      s          ds\n\n');






%% Start iteration with the secant method
%
while(abs(ds) > epsi) & (it < itmax)
   it = it + 1;
   y0 = [0.0 0.0 s1];
   [x,y] = ode45(@blasius,xspan,y0,options);
   fi1 = y(end,2) - 1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s0 = s1;
   s1 = s1 + ds;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s1,ds);
end




%% Compute a table for f, f' and f''
xspan = (x0:0.25:x1);
y0 = [0.0 0.0 s1];
[x,y] = ode45(@blasius,xspan,y0,options);
fprintf('\n         eta        f          f''        f"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[x y]');
% Plotting f' and f" as a function of eta
xspan = [x0 x1];
y0 = [0.0 0.0 s1];
[x,y] = ode45(@blasius,xspan,y0,options);

h=plot(y(:,2),x, y(:,3),x,'-.');


%% Legends, labels and title
grid on
xlabel('f'' , f"');
ylabel('\eta');
h3=legend('f''','f"');

set(h3,'box','off');
title('Solution of Blasius equation');
