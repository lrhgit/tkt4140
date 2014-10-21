%================= Program nlsor =========================
% Using SOR for solving the following ODE:
%
%          y''(x) = 1.5*y(x)^2
%           y(0) = 4, y(1) = 1
%
% The scheme is given in eq. 3.7.8 in the compendium
% w is the relaxation factor and h the steplenght 
% Initial values : y = 0 in the whole interval except for
% the boundary values y(0) = 4 and y(1 = 1
%=========================================================
clear;
h = input(' Steplenght h = ?');
w = input('Relaxation factor w = ? ');
ni = 1/h ; % No. of intervals
% Choose h so that ni is an integer
np = ni + 1 ; % No. of points
fprintf('\n  Steplength h = %9.3e \n',h);
fprintf('  Relaxation factor w  = %6.3f \n',w);
y = zeros(np,1); % Initial values
y(1) = 4.0; y(np) = 1.0;
epsi = 1.0e-5; it = 0; test = 1.0; h2 = h*h;
while (test > epsi)
   abserr = 0.0;
   for k = 2 : ni 
      f = y(k-1)- y(k)*(2 + 1.5*h2*y(k)) + y(k+1);
      df = - (2 + 3*h2*y(k));
      dy = - w*f/df;    
      y(k) = y(k) + dy;      
      abserr = max(abserr, abs(dy));
   end    
   test = abserr;   
   it = it + 1;
end
fprintf('\n No. of iterations = %8.0f \n',it); 
