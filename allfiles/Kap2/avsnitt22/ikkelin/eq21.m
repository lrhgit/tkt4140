% ========= Program eq21 =========
%  Finner s = y'(0) for eksemplet
%  i avsnitt 2.2 i kompendiet
%  
%  Ligningen er gitt ved:
%  y''(x) = a*y^2 , a > 0
%  y(0) = 6/a, y(1) = 3/(2*a)
%  
%  I programmet setter vi y(1)= y, y(2) = y'
% 
% s = y'(0) = y(2)(0)
% Bruker funksjonen fcn21.m
%
% Skriv inn a samt startverdier s0 og s1
%
clear; clear global a;
global a;
a = 10.0;
xspan = [0 1];
y0 = 6/a;
y1 = 3/(2*a);
% Velg startverdier
s0 = -8.0;
s1 = -5.0;
% Beregner fi0
ys = [y0 s0];
options = odeset('RelTol',1.0e-7);
[x,y] = ode45(@fcn21,xspan,ys,options);
fi0 = y(end,1) - y1 ;
% Startverdier for iterasjonen
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
% Skriv overskrift for tabell
fprintf('        itr.      s          ds\n\n');
% Starter iterasjon
while(abs(ds) > epsi) && (it < itmax)
   it = it + 1;
   ys = [y0 s1];
   [x,y] = ode45(@fcn21,xspan,ys,options);
   fi1 = y(end,1) - y1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
% Beregner en tabell for y og y' 
% xspan = (0:0.1:1);
% y0 = [4.0 s];
% [x,y] = ode45('fcn21',xspan,y0,options);
% fprintf('\n         x       y          y'' \n\n');
% fprintf(' %12.2f %10.6f %10.6f \n',[x y]');
% % Plotter y 
% xspan = [0 1];
% y0 = [4.0 s];
% [x,y] = ode45('fcn21',xspan,y0,options);
% plot(x,y(:,1))
% grid on
% xlabel('x','FontSize',14,'FontWeight','Bold')
% ylabel('y','Fontsize',14)
% title('Eksempel i avsnitt 2.2','Fontsize',14)
% %legend('y','y''')