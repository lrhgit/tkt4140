function fsksec2 
% Solution of the Falkner-Skan equation using a
% shooting technique. (See example 2.5 in the compendium)
% The equation is given by:
%    f'''(eta) + f(eta)*f''(eta)+ beta*(1- (f')^2) = 0
%    f(0) = 0, f'(0) = 0, f'(etainf)=1
% In the program we put y(1)= f, y(2) = f'
% and y(3) = f'' using x instead of eta
% The program reads the beta, s0, s1 and etainf. 
%  s = f''(0) = y(3)(0)
%
% Using function fcnfsk
%
clear; clear global beta;
global beta;
beta = input(' beta = ?');
s0 = input(' s0 = ?');
s1 = input(' s1 = ?');
x1 = input(' etainf = ?');
% === Print beta, s0 ,s1 and etainf
fprintf('        beta = %8.4f\n',beta);
fprintf('        s0 = %8.5f s1 = %8.5f  etainf = %7.3f \n\n',s0,s1,x1);
x0 = 0; 
xspan = [x0 x1];
% Compute fi0
y0 = [0.0 0.0  s0];
[x,y] = ode45(@fcnfsk,xspan,y0);
fi0 = y(end,2) -1;
% Initial values for the iteration
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
% Heading of table
fprintf('        itr.      s          ds\n\n');
% Start iteration
while(abs(ds) > epsi) & (it < itmax)
   it = it + 1;
   y0 = [0.0 0.0 s1];
      [x,y] = ode45(@fcnfsk,xspan,y0,options);
   fi1 = y(end,2) -1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
% Compute a table of f, f' og f'' 
dx = 0.25;
xspan = (x0:dx:x1);
if(mod(x1,dx)~=0)
   xspan = [xspan x1];
end
y0 = [0.0 0.0 s];
[x,y] = ode45(@fcnfsk,xspan,y0,options);
fprintf('\n         eta        f          f''        f"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[x y]');
% Plotting f' and f" as a function of eta
xspan = [x0 x1];
y0 = [0.0 0.0 s];
[x,y] = ode45(@fcnfsk,xspan,y0,options);
plot(y(:,2),x)
hold on
plot(y(:,3),x,'-.')
grid on
ylabel('\eta','FontSize',14,'FontWeight','Bold','Rotation',0)
xlabel('f'' , f"','Fontsize',14)
title('Solution of the Falkner-Skan equation','Fontsize',14)
legend('f''','f"')
hold off
shg
%-----------------------------------------------
function yd = fcnfsk(x,y)
% Used by fsksec
global beta;
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -(y(1)*y(3) + beta*(1 - y(2)^2));