function [x,y]=falknerSkanPenta(beta,h,xmax) 
%
% falkerSkanPenta computes the solution of the Falkner-Skan
% equation using finite differences:
%            f''' + f*f'' + beta*(1 - (f')^2) = 0
%         with boundary conditions :
%                f(0) = f'(0) = 0
%          f'(eta) -> 1.0, eta  -> infinity
% This version have a penta-diagonal matrix of coefficients.
% Input-variables : h and iprt where h is the grid-
% spacing and iprt is an output-control parameter.
% iprt = n , n > 0, prints results for every n gridpoint.
% iprt = 0 prints output for only the first and last gridpoint
%
% Uses function penta and svalue
% 
% betasep = -0.19883768;
% beta = input(' beta = ?');
% while(beta < betasep) | (beta > 1.999)
%    fprintf('beta = %7.3e is invalid!\n',beta);
%    disp(' Try again !');
%    beta = input(' beta = ?');
% end

%etainf = svalue(beta);
%h = input(' Input grid-spacing h : ');
% iprt = input(' Input print-control iprt :');
iprt=1;
iprt = abs(iprt);
ni = round(xmax/h); % Number of intervals
etainf = ni*h; % Recompute etainf (safe play)
n = ni - 1 ; % Number of equations

x = linspace(h,xmax,ni)';

itmax = 10 ; epsi = 1.0e-5;


% --- Initialize
y = zeros(n,1);  
e = -ones(n,1); f = -e; % Not destroyed in Penta
b = zeros(n,1); a = b; c = b ; d = b; ytmp = b;


% === Initialize y using Pohlhausen polynomials
t1 = h/etainf;
for k = 1 : n
   ksi = t1*k;
   ksi2 = ksi*ksi;
   y1 = 1.0 + ksi2*(ksi*0.2 -0.5);
   y2 = beta*(0.5 - ksi + ksi2*(0.75 - 0.2*ksi))*etainf^2/6;
   y(k) = etainf*ksi2*(y1 + y2);
end

% if iprt > 0
%    peta = iprt*h;
%    fprintf('   Print-spacing  ............................%10.3e\n',peta);
% else
%    fprintf('   Print-spacing  ............................%10.3e\n',etainf);
% end
% fprintf('\n           Number of intervals = %10.0f \n',ni);
% fprintf('           Max. value of eta   = %10.4f \n',etainf);
% fprintf('\n      Iteration \n');
% fprintf('            no.     mean error \n');

% -------------------------------
%     Start of iteration loop 
% -------------------------------

it = 0; merror = 1.0;
bh = beta*h; bh3 = beta*h^3;
tic
while (it <= itmax)& (merror > epsi)
   it = it + 1;
   for k = 2 : n-1
      a(k) = 2*(h*y(k) + 1) + bh*(y(k+1) - y(k-1));
      b(k) = 2*h*(y(k+1) - 4*y(k) + y(k-1));
      c(k) = 2*(h*y(k) - 1) - bh*(y(k+1) - y(k-1));
      d(k) = 2*h*y(k)*(y(k+1) - 2*y(k) + y(k-1))...
             -2*bh3 - 0.5*bh*(y(k+1) - y(k-1))^2;
   end
   a(n) = 2*(h*y(n) + 1) + bh*(y(n) - y(n-1) + h);
   b(1) = 2*h*(y(2) - 4*y(1)) - 1;
   c(1) = 2*(h*y(1)-1) - bh*y(2);
   b(n) = 2*h*(y(n-1) - 2*y(n)) + 2*h^2 - 1 -bh*(y(n)-y(n-1) + h); 
   c(n-1) = 2*h*y(n-1)-1 - bh*(y(n) - y(n-2));
   d(1) =   2*h*y(1)*(y(2) - 2*y(1))- 2*bh3 - 0.5*bh*y(2)^2;
   d(n-1) = - h + 2*h*y(n-1)*(y(n) - 2*y(n-1) + y(n-2)) ...
            -2*bh3 -0.5*bh*(y(n) - y(n-2))^2;
   d(n) = 2*h*y(n)*(y(n-1) - y(n))-1.5*bh3 - 0.5*bh*(y(n) - y(n-1))^2;
   
   ytmp = penta(e,a,b,c,f,d); % Solve system of equations
   merror = sum(abs(y-ytmp))/n; %Compute mean error
   fprintf('     %10.0f   %12.3e  \n',it, merror);
   y = ytmp; % Update the y-vector
end
if it >= itmax
   fprintf('\n  *** Max. number of iterations! ***\n');  
end;
toc


%% -----------------------
%   Printing of tables 
%-----------------------
yend = y(end)+ h;
y =[0;y;yend];
np = length(y);
for k = 2:ni
   yd(k) = (y(k+1) - y(k-1))*0.5/h;
end
yd(1) = 0; yd(np) = 1.0;
for k = 2:ni
   y2d(k) = (yd(k+1) - yd(k-1))*0.5/h;
end
y2d(1) = 2*y(2)/h^2;
y2d(np) = (1- yd(np-1))*0.5/h;

y =[yd;y2d];


%%-------------------------------------------------
function etamax = svalue(beta)
% Estimate of etamax for a given value of beta
% See appendix 3,part 3.
p1 = [0.633 -1.68 5.76] ; p2 = [36.76  2.0 5.87];
p3 = [0.125  -0.9 5.463]; 
if ( beta <= 1.0)
   if ( beta >= 0.0 ) % 0 <= beta <= 1.0
      p = polyval(p1,beta);
      etamax = round(10.0*p)*0.1;
   else
      p = polyval(p2,beta); % betasep <= beta < 0
      etamax = round(10.0*p)*0.1;
   end
else
   p = polyval(p3,beta); % 1 < beta <= 1.999
   etamax = round(10.0*p)*0.1;
end

%--------------------------------------------------
function x = penta(e,a,b,c,f,d)     
n = length(b);
x = zeros(size(b));
% === Elimination
for k = 2 : n-1 
   q = a(k)/b(k-1);
   b(k) = b(k) - q*c(k-1);
   c(k) = c(k) - q*f(k-1);
   d(k) = d(k) - q*d(k-1);
   q = e(k+1)/b(k-1);
   a(k+1) = a(k+1) - q*c(k-1);
   b(k+1) = b(k+1) - q*f(k-1);
   d(k+1) = d(k+1) - q*d(k-1);
end
% === Backsubstitution
q = a(n)/b(n-1);
b(n) = b(n) - q*c(n-1);
x(n) = (d(n) - q*d(n-1))/b(n);
x(n-1) = (d(n-1) - c(n-1)*x(n))/b(n-1);
for k = n-2 : -1 :1   
   x(k) = (d(k) - f(k)*x(k+2) - c(k)*x(k+1))/b(k); 
end





