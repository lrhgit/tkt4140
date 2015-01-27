%===================== Wedge =========================
% Programmet løser ligningen for varmeledning i en 
% kile sammensatt av to materialer med bruk av skyteteknikk.
% Ligningen er gitt ved:
%     (x*theta'(x)/beta^2)' = theta(x)
%     beta^2 = 2 for x<= 0.5, beta^2 = 8 for x > 0.5
%     theta'(0) = 2*theta(0), theta(1) = 1
% Løser systemet med bruk av RK4C
%========================================================
clear;
dx = 0.01;
n = 1.0/dx; th = zeros(n+1,1); dth = th;
y = zeros(2,1);
% ==== Tipper initialverdier s0 og s1:
s0 = input('s0 = ?');
s1 = input('s1 = ?');
fprintf('  Initialverdier: s0 = %7.3f  s1 = %7.3f \n',s0,s1);

% === Beregner fi0 ===
y(1) = s0 ; y(2) = 0;
for k = 1: n
   x = (k-1)*dx;  
   y = RK4C('fcnwed',x,y,dx);
end
fi0 = y(1) - 1;
fprintf('  fio = %10.3e \n',fi0);
% === Beregner fi1 ===
 y(1) = s1 ; y(2) = 0 ;
for k = 1: n 
   x = (k-1)*dx;  
   y = RK4C('fcnwed',x,y,dx);
end
fi1 = y(1) - 1;
fprintf('  fi1 = %10.3e \n',fi1);
% === Beregner s* ===
sstar = (fi1*s0 - s1*fi0)/(fi1 - fi0) ;
fprintf('  s*  = %10.3e \n',sstar);
% === Endelig beregning med s* ===
y(1) = sstar ; y(2) = 0 ;
for k = 1: n 
   x = (k-1)*dx;  
   y = RK4C('fcnwed',x,y,dx);
   th(k+1) = y(1);
   dth(k+1) = y(2);
end
th(1) = sstar;
dth(1) = 2*th(1);
b2 = 2.0;
for k = 2 : n+1
   x = (k - 1)*dx;
   if ( x > 0.5 )
      b2 = 8.0;
   end
   dth(k) = b2*dth(k)/x;
end

x = [0 : dx : 1.0]';
% === Tabell ===
fprintf('\n    x           theta       dtheta \n\n'); 
fprintf(' %7.3f   %10.4e  %10.4e \n',[x th dth]'); 


