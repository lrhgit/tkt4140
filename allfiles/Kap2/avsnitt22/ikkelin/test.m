% === Program test ===
% Programmet løser randverdi-problemet i 
% avsnitt 2.2 i kompendiet ved bruk av skyteteknikk.
% Bruker sekantmetoden for nullpunktbestemmelse.
% Ligning : y''(x)= 6*y(x)^2
%           y(0) = 1, y(1) = 0.25
% 
clear; 
xspan = [0.0 1.0];
% ----Tipper to verdier s0 og s1 for y'(x0)
s0 = -8.0; s1 = -9.0;
% ---- Beregner fi0
y0 = [1.0; s0]; % Startverdier
[x,y] = ode45(@fcn22,xspan,y0);
fi0 = y(end,1) -0.25;
% ----- Startverdier for iterasjonen
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
fprintf('        itr.      s          ds\n\n');
% ===== Starter iterasjon =====
while(abs(ds) > epsi) & (it < itmax)
   it = it + 1;
   y0 = [1.0; s1];
   [x,y] = ode45(@fcntest,xspan,y0,options);
   fi1 = y(end,1) -0.25;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s0 = s1;
   s1 = s1 + ds;
   fi0 = fi1;
   if abs(s1) > 1 
      ds = ds/s1;
   end;
   fprintf('%10d %12.6f %12.3e\n',it,s1,ds);
end
% ---- Beregner en tabell for y og y'
xspan = [0:0.025:1.0];
y0 = [1.0; s1];
options = odeset('RelTol',1.0e-5);
[x,y] = ode45(@fcntest,xspan,y0,options);
fprintf('\n           x         y           y''\n\n');
fprintf(' %12.3f %12.6f %12.6f\n',[x y]');
