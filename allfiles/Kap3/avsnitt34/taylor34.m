<<<<<<< HEAD
%=================== Program taylor34  ==============
% Solves the example in section 3.4 where the
% differential equation is given in eq. 4.1, chap. 3
% The difference eqation is linearized using 
% Taylor-expansionng
% Solution by tdma
%====================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.05; % steplength
ni = 1/h; % No. of intervals
% Choose h so that ni is an integer.
n = ni-1; % No. of equations
fac = 3.0*h*h;
nitr = 6; % No. of iterations
a = ones(n,1) ; % subdiagonal
% a  is no destroyed in the elimination 
%--- Initialize
c = zeros(n,1); % superdiagonal 
ym = c; b = ym; d = ym; 



fprintf('        Itr.      max. deviat.  \n');

it = 0; itmax = 15; dymax = 1.0; RelTol = 1.0e-10;

while (dymax > RelTol) & (it < itmax)
   it = it +1;	
   c = ones(n,1);
   b = -(2.0 + fac*ym); % main diagonal
   d = -(fac*0.5)*ym.^2 ; % right hand side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   ym1 = tdma(a,b,c,d); % Solve the system
   dymax = max(abs((ym1-ym)./ym1));% Compute rel. deviation
   ym = ym1; % Update the y-values
   fprintf(' %10d     %12.3e \n',it,dymax);
end

%---- Print y and relative error ----
xv = (h:h:1.0-h)';
ya = 4.0./(1 + xv).^2; % Analytical solutionk lsning
feil = abs((ym1 - ya)./ya);
fprintf('\n      x         y     rel. error  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv ym1 feil]');

=======
%=================== Program taylor34  ==============
% Solves the example in section 3.4 where the
% differential equation is given in eq. 4.1, chap. 3
% The difference eqation is linearized using 
% Taylor-expansionng
% Solution by tdma
%====================================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

h = 0.05 ; % steplength
ni = 1/h; % No. of intervals
% Choose h so that ni is an integer.
n = ni-1; % No. of equations
fac = 3.0*h*h;
nitr = 6; % No. of iterations
a = ones(n,1) ; % subdiagonal
% a  is no destroyed in the elimination 
%--- Initialize
ym = zeros(n,1); b = ym; d = ym; 
c = ones(n,1);


fprintf('        Itr.      max. deviat.  \n');

it = 0; itmax = 15; dymax = 1.0; RelTol = 1.0e-10;

while (dymax > RelTol) & (it < itmax)
   it = it +1;	
   b = -(2.0 + fac*ym); % main diagonal
   d = -(fac*0.5)*ym.^2 ; % right hand side 
   d(n) = d(n)- 1.0;
   d(1) = d(1) - 4.0;
   ym1 = tdma(a,b,c,d); % Solve the system
   dymax = max(abs((ym1-ym)./ym1));% Compute rel. deviation
   ym = ym1; % Update the y-values
   fprintf(' %10d     %12.3e \n',it,dymax);
end

%---- Print y and relative error ----
xv = (h:h:1.0-h)';
ya = 4.0./(1 + xv).^2; % Analytical solutionk lsning
feil = abs((ym1 - ya)./ya);
fprintf('\n      x         y     rel. error  \n\n')
fprintf('  %7.3f %10.5f  %10.2e \n',[xv ym1 feil]');

>>>>>>> master
plot(xv,ym1,xv,ya)