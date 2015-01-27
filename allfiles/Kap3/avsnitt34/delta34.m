%====================== Program delta34  =====================
% Solves the example in section 3.4 where the
% differential equation is given in equation 4.1,chap. 3
% The difference equation is linearized using Taylor-expansion
% The equations are in delta-form.
% Solution by tdma.
%=============================================================
clear
h = 0.05; % steplenght
ni = 1/h; % No. of intervals
% choose h so that ni is an integer
n = ni - 1; % No. of equations
fac = 3.0*h*h;
%--- Initialize
a = ones(n,1) ; % subdiagonal
c = a; % superdiagonal
b = a; d = a; dy = a;
y = zeros(n,1); 
fprintf('        Itr.      max. deviat.  \n');
it = 0; itmax = 10; dymax = 1.0; RelTol = 1.0e-5;
while (dymax > RelTol) & (it < itmax)    
  it = it + 1;
   b = -(2.0 + fac*y); % main diagonal
   d = (fac*0.5)*y.^2 ; % right hand side 
   for j = 2:n-1
      d(j) = d(j) -(y(j+1)-2*y(j) + y(j-1));
   end
   d(n) = d(n)- (1.0- 2*y(n) + y(n-1));
   d(1) = d(1) - (y(2)-2*y(1) + 4.0);
   dy = tdma(a,b,c,d); % Solve the system
   y = y + dy; % Update the y-values
   dymax = max(abs((dy)./y));% Compute relative deviation 
   fprintf(' %10d     %12.3e \n',it,dymax);
end
%---- Ptint y and relative error ----
xv = (h:h:1.0-h)';
ya = 4.0./(1 + xv).^2; % Analytical solution
feil = abs((y - ya)./ya);
fprintf('\n      x         y      rel. error  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv y feil]');
