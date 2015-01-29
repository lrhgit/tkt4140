% ===================== Program bladoug ===========================
% This program computes the solution of the Blasius
% equation using finite differences and the Douglas algorithm:
%            f''' + f*f'' = 0
%         with boundary conditions :
%                f(0) = f'(0) = 0
%          f'(eta) -> 1.0, eta  -> infinity
%
% The system of equations is solved using function douglas which is
% a stripped versjon of function bitri.
%
% Input-variables : deta and iprt where deta is the grid-
% spacing and iprt is an output-control parameter.
% iprt = n , n > 0, prints results for every n gridpoint
%                   including the last point.    
% iprt = 0 prints output for only the first and last gridpoint
% =================================================================
clear
etainf = 5.8;
deta = input(' Input grid-spacing deta : ');
iprt = input(' Input print-control iprt :');
iprt = abs(iprt);
nval = round(etainf/deta); % Number of intervals
etainf = nval*deta;
neq = nval - 1 ; % Number of equations
% ---- Allocation and initial values----
c2 = ones(neq,1); a2 = -c2; b3 = -2*c2; b1 = -2*deta*c2;
d1 = zeros(neq,1); a3 = d1;b2 = d1; b4 = d1; c3 = d1; d2 = d1;
f = zeros(neq,1); g = f; u = f; 
%
itmax = 10; it = 0; merror = 1; epsi = 1.0e-5;
disp('        ****************************************');
disp('        *         BLASIUS EQUATION             *');
disp('        *                                      *');
disp('        *         Douglas algorithm            *');
disp('        ****************************************');
fprintf('\n   Max. number of iterations, itmax ..........%2.0f\n',itmax);
fprintf('   Max. allowable error in iteration, epsi ...%10.3e\n',epsi);
fprintf('   Grid-spacing, deta ........................%10.3e\n',deta);

% === Initialize f and f'= g using Pohlhausen polynomials
t1 = deta/etainf;
for k = 1 : neq
   ksi = t1*k; ksi2 = ksi*ksi;
   f(k) = etainf*ksi2*(1.0 +ksi2*(ksi*0.2 -0.5));
   g(k) = ksi*(2.0 + ksi2*(ksi - 2.0));
end
if iprt > 0
   peta = iprt*deta;
   fprintf('   Print-spacing  ............................%10.3e\n',peta);
else
   fprintf('   Print-spacing  ............................%10.3e\n',etainf);
end
fprintf('\n           Number of grid points = %10.0f \n',nval);
fprintf('           Max. value of eta     = %10.4f \n',etainf);
disp (' ');
disp('      Iteration');
disp('            no.   mean err.');
disp(' '); 
fac = deta*0.5; 
itmax = 10; it = 0; merror = 1; epsi = 1.0e-5;
while (merror > epsi) & (it <= itmax)% Start of iteration loop
   it = it + 1;
   for k = 2 : neq-1;
       a3(k) = 1-fac*f(k);
       b4(k) = fac*(g(k+1) - g(k-1));
       c3(k) = 1 + fac*f(k);
       d2(k) = fac*f(k)*(g(k+1) - g(k-1));
   end
   b4(1) = fac*g(2);
   c3(1) = 1 + fac*f(1);
   d2(1) = fac*f(1)*g(2);
   a3(neq) = 1 -fac*f(neq);
   b2(neq) = 1;
   b4(neq) = fac*(1 - g(neq-1));
   d1(neq) = -deta;
   d2(neq) = -(1 + fac*f(neq)*g(neq-1));
   % Solve system of equations
   [u,f] = douglas(a2,a3,b1,b2,b3,b4,c2,c3,d1,d2); 
   merror = sum(abs(u-g))/neq;
   g = u;
   fprintf('     %10.0f  %10.3e \n',it,merror);
end % end of iteration loop
if it > itmax
    disp(' *** Max. number of iterations! ***');
end; 
% -----------------------
%    Printing of tables 
% -----------------------
fprintf('\n  eta       f             f''            f" \n\n');
fetainf = f(neq) + deta;
f = [0; f;fetainf]; g = [0;g;1];
gd = zeros(length(g),1);
for k = 2:neq
    gd(k) = (g(k+1) - g(k-1))*0.5/deta;
end
gd(1) = 2*f(2)/deta^2;
eta = (0 : deta : etainf);
npr = length(eta);
if iprt > 0
    for k = 1 : iprt: npr
        fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(k),f(k),g(k),gd(k));
    end
    if k ~= npr
        k = npr;
        fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(k),f(k),g(k),gd(k));
    end
else
    k = 1;
    fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(k),f(k),g(k),gd(k));
    k = npr;
    fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(k),f(k),g(k),gd(k));
end





