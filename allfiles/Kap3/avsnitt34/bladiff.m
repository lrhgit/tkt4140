% ================ Program bladif ===========================
%
% This program computes the solution of the Blasius
% equation using finite differences:
%            f''' + f*f'' = 0
%         with boundary conditions :
%                f(0) = f'(0) = 0
%          f'(eta) -> 1.0, eta  -> infinity
% Input-variables : deta and iprt where deta is the grid-
% spacing and iprt is an output-control parameter.
% iprt = n , n > 0, prints results for every n gridpoint.
% iprt = 0 prints output for only the first and last gridpoint
%
% ============================================================
clear
etainf = 5.8;
deta = input(' Input grid-spacing deta : ');
iprt = input(' Input print-control iprt :');
iprt = abs(iprt);
np = round(etainf/deta); % Number of gridpoints
etainf = np*deta;
neq = np - 1 ; % Number of equations
f = zeros(np,1); g = f; gd = f;
b = zeros(neq,1); a = b; c = b ; d = b; x = b;
itmax = 15; epsi = 1.0e-5;
disp('        ****************************************');
disp('        *         BLASIUS EQUATION             *');
disp('        *                                      *');
disp('        *         Difference method            *');
disp('        ****************************************');
fprintf('\n   Max. number of iterations, itmax ..........%2.0f\n',itmax);
fprintf('   Max. allowable error in iteration, epsi ...%10.3e\n',epsi);
fprintf('   Grid-spacing, deta ........................%10.3e\n',deta);

% === Initialize f and f' using Pohlhausen polynomials
t1 = deta/etainf;
for k = 1 : np
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
f2wold = 2.0*f(1)/deta^2;
fprintf('\n           Number of grid points = %10.0f \n',np);
fprintf('           Max. value of eta     = %10.4f \n',etainf);
fprintf('           Estimate of f"(0)     = %10.4f \n\n',f2wold);
disp (' ');
disp('      Iteration');
disp('            no.   s = f"(0)          ds ');
disp(' ');

% -------------------------------
%     Start of iteration loop 
% -------------------------------
it = 0;  fac = deta*0.5; df2 = 1.0;
while (df2 > epsi) & (it <= itmax)
   it = it + 1;
   a = 1.0 - fac*f;
   b = -2.0*ones(neq,1);
   c = 1.0 + fac*f;
   d(neq) = -c(neq);
   x = tdma(a,b,c,d); % Solve system of equations
   for k = 1 : neq
      g(k) = x(k);
   end
   % === Compute new f by the trapezoidal method ===
   f(1) = fac*g(1);
   s = f(1);
   for k = 2 : np 
      s = s + fac*(g(k) + g(k-1)); 
      f(k) = s;  
   end
   f2new = 2*f(1)/deta^2;
   df2 = abs(f2new - f2wold);
   fprintf('     %10.0f  %10.4e    %10.3e \n',it,f2new,df2);
   f2wold = f2new;
end
if it > itmax
    disp(' *** Max. number of iterations! ***');
end;
t2 = 0.5/deta;
for k = 2 : neq
   gd(k) = (g(k+1) - g(k-1))*t2;
end
gd(1) = g(2)*t2;
gd(np) = (1.0 - g(neq))*t2;
% -----------------------
%    Printing of tables 
% -----------------------
fprintf('\n  eta       f             f''            f" \n\n');
f = [0; f]; g = [0;g]; gd = [f2new;gd];
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
    fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(1),f(1),g(1),gd(1));
    fprintf('%6.3f  %10.5e  %10.5e  %10.5e \n',eta(npr),f(npr),g(npr),gd(npr));
end





