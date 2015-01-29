%===================== wedge45 =========================
% Programmet løser ligningen for varmeledning i en 
% kile sammensatt av to materialer med bruk av skyteteknikk.
% Ligningen er gitt ved:
%     (x*theta'(x)/beta^2)' = theta(x)
%     beta^2 = 2 for x<= 0.5, beta^2 = 8 for x > 0.5
%     theta'(0) = 2*theta(0), theta(1) = 1
% Løser systemet med bruk av ode45
%========================================================
clear;
% ==== Tipper initialverdier s0 og s1:
s0 = 0.06;
s1 = 0.1;
fprintf('  Initialverdier: s0 = %7.3f  s1 = %7.3f \n',s0,s1);
xint = [0 1];
options = odeset('RelTol',1.0e-6);
% === Beregner fi0 ===
y0 = [s0 ; 0];
[x,y] = ode45(@fcnwed,xint,y0,options);
fi0 = y(end,1) - 1;
fprintf('  fio = %10.3e \n',fi0);

% === Beregner fi1 ===
y0 = [s1 ; 0];
[x,y] = ode45(@fcnwed,xint,y0,options);
fi1 = y(end,1) - 1;
fprintf('  fi1 = %10.3e \n',fi1);

% === Beregner s* ===
sstar = (fi1*s0 - s1*fi0)/(fi1 - fi0) ;
fprintf('  s*  = %10.3e \n',sstar);

% === Endelig beregning med s* ===
y0 = [sstar ; 0];
xint = (0: 0.05 :1)';
[x,y] = ode45(@fcnwed,xint,y0,options);
n = length(y);
dth = zeros(n,1);
dth(1) = 2*y(1,1);
b2 = 2.0;
for k = 2 : n
    xk = x(k);
   if ( xk > 0.5 )
      b2 = 8.0;
   end
   dth(k) = b2*y(k,2)/xk;
end
% === Tabell ===
fprintf('\n    x           theta       dtheta \n\n'); 
fprintf(' %7.3f   %10.5e   %10.5e \n',[x y(:,1) dth]'); 


