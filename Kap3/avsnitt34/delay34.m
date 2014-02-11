%================ Program delay34  ===============
% Solves the example in section 3.4 where the 
% differential equation is equation 4.1 in chap.3
% The difference equation is using the method of 
% delayed coefficients.
% The system is is solved by function tdma
%=================================================
clear
h = 0.05; % steplenght
n = round(1.0/h)-1; % No. of equations
fac = (3.0/2.0)*h*h;
nitr = 12; % No. of iterations
%--- a  is not destroyed in the elimination
a = ones(n,1) ; % subdiagonal
%--- Initializing
c = zeros(n,1); % superdiagonal 
ym = c; b = ym; d = ym; 
fprintf('        Itr.      max. deviat.  \n');
for k = 1:nitr 
   c = ones(n,1);
   b = -(2.0 + fac*ym); % main diagonal
   d(n) = - 1.0;
   d(1) =  - 4.0;
   ym1 = tdma(a,b,c,d); % Solve the system
   dymax = max(abs((ym1-ym)./ym1));% Compute relative deviation
   ym = ym1;
   fprintf(' %10d     %12.3e \n',k,dymax);
 end
 %---- Print  y and error ----
 xv = (h:h:1.0-h)';
 ya = 4.0./(1 + xv).^2; % Analytical solution
 feil = abs((ym1 - ya)./ya);
fprintf('\n      x         y      rel. error  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv ym1 feil]');
