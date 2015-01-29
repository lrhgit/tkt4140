% ========================= Blapenta ==========================
%
% Blapenta computes the solution of the Blasius
% equation using finite differences:
%            f''' + f*f'' = 0
%         with boundary conditions :
%                f(0) = f'(0) = 0
%          f'(eta) -> 1.0, eta  -> infinity
% This version have a penta-diagonal matrix of coefficients.
% Input-variables : h and iprt where h is the grid-
% spacing and iprt is an output-control parameter.
% iprt = n , n > 0, prints results for every n gridpoint.
% iprt = 0 prints output for only the first and last gridpoint
%
% =============================================================
clear
etainf = 5.8;
h = input(' Input grid-spacing h : ');
iprt = input(' Input print-control iprt :');
iprt = abs(iprt);
ni = round(etainf/h); % Number of intervals
etainf = ni*h; % Recompute etainf (safe play)
neq = ni - 1 ; % Number of equations
% --- Initialize
y = zeros(neq,1);  
e = -ones(neq,1); f = -e; % Not destroyed in Penta
b = zeros(neq,1); a = b; c = b ; d = b; x = b;
itmax = 10 ; epsi = 1.0e-5;
fprintf('        ****************************************\n');
fprintf('        *           BLASIUS EQUATION           *\n');
fprintf('        *                                      *\n');
fprintf('        *         Difference method with       *\n');
fprintf('        *         pentadiagonal matrix         *\n');
fprintf('        *                                      *\n');
fprintf('        ****************************************\n');
fprintf('\n   Max. number of iterations, itmax ..........%2.0f\n',itmax);
fprintf('   Max. allowable error in iteration, epsi ...%10.3e\n',epsi);
fprintf('   Grid-spacing, deta ........................%10.3e\n',h);

% === Initialize y using Pohlhausen polynomials
t1 = h/etainf;
for k = 1 : neq
   ksi = t1*k; ksi2 = ksi*ksi;
   y(k) = etainf*ksi2*(1.0 + ksi2*(ksi*0.2 -0.5));
end
if iprt > 0
   peta = iprt*h;
   fprintf('   Print-spacing  ............................%10.3e\n',peta);
else
   fprintf('   Print-spacing  ............................%10.3e\n',etainf);
end
fprintf('\n           Number of intervals = %10.0f \n',ni);
fprintf('           Max. value of eta   = %10.4f \n',etainf);
fprintf('\n      Iteration \n');
fprintf('            no.     mean error \n');

% -------------------------------
%     Start of iteration loop 
% -------------------------------
it = 0; merror = 1.0;
while (it <= itmax)& (merror > epsi)
   it = it + 1;
   a = 2*(h*y + 1);
   c = 2*(h*y - 1);
   for k = 2 : neq-1
      b(k) = 2*h*(y(k+1) - 4*y(k) + y(k-1));
      d(k) = 2*h*y(k)*(y(k+1) - 2*y(k) + y(k-1));
   end
   b(1) = 2*h*(y(2) - 4*y(1)) - 1;
   d(1) = 2*h*y(1)*(y(2) - 2*y(1));
   b(neq) = 2*h*(y(neq-1) - 2*y(neq)) + 2*h^2 - 1;
   c(neq-1) = 2*h*y(neq-1)-1;
   d(neq-1) = -h + 2*h*y(neq-1)*(y(neq) - 2*y(neq-1) + y(neq-2));
   d(neq) = 2*h*y(neq)*(y(neq-1) - y(neq));
   
   x = penta(e,a,b,c,f,d); % Solve system of equations
      
   % --- Compute mean error
   merror = sum(abs(y-x))/neq;
   fprintf('     %10.0f   %12.3e  \n',it, merror);
   
   y  = x ; % Update the y-vector
end
if it >= itmax
   fprintf('\n  *** Max. number of iterations! ***\n');  
end;
% -----------------------
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
fprintf('\n  point   eta        f              f''             f" \n\n');
k = 1;
eta = 0;
s1 = sprintf(' %5.0f  %6.3f ',k-1,eta);
s2 = sprintf(' %13.5e  %13.5e  %13.5e \n',y(1),yd(k),y2d(k));
fprintf([s1 s2]);
if iprt > 0
   while k<= np - iprt
      k = k + iprt;
      eta = (k-1)*h;
      s1 = sprintf(' %5.0f  %6.3f ',k-1,eta);
      s2 = sprintf(' %13.5e  %13.5e  %13.5e \n',y(k),yd(k),y2d(k));
      fprintf([s1 s2]);
    end
end
if k ~= np
   k = np;
   eta = (k-1)*h;
   s1 = sprintf(' %5.0f  %6.3f ',k-1,eta);
   s2 = sprintf(' %13.5e  %13.5e  %13.5e \n',y(k),yd(k),y2d(k));
   fprintf([s1 s2]);
end





