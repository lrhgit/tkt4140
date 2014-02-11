function [x,y]=falknerSkanSecant(beta)
% Solution of the Falkner-Skan equation using a
% shooting technique. (See example 2.5 in the compendium)
% The equation is given by:
%    f'''(eta) + f(eta)*f''(eta)+ beta*(1- (f')^2) = 0
%    f(0) = 0, f'(0) = 0, f'(etainf)=1
% In the program we put y(1)= f, y(2) = f'
% and y(3) = f'' using x instead of eta
% The program reads the beta-value and uses the 
% function svalue to compute an estimate of s0,s1 and etainf
%  s = f''(0) = y(3)(0)
%
% Using functions fcnfsk and svalue
%
global beta

% === Compute initial values
[s0,s1,x1] = svalue(beta); % x1 = etainf
% === Print beta, s0 ,s1 and etainf
% fprintf('        beta = %8.4f\n',beta);
% fprintf('        s0 = %8.5f s1 = %8.5f  etainf = %7.3f \n\n',s0,s1,x1);
x0 = 0; 
xspan = [x0 x1];
% Compute fi0
y0 = [0.0 0.0  s0];
[x,y] = ode45(@fcnfsk,xspan,y0);
fi0 = y(end,2) -1;
% Initial values for the iteration
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-7);
% Heading of table
fprintf('        itr.      s          ds\n\n');

%% Start iteration
tic
while(abs(ds) > epsi) & (it < itmax)
   it = it + 1;
   y0 = [0.0 0.0 s1];
      [x,y] = ode45(@fcnfsk,xspan,y0,options);
   fi1 = y(end,2) -1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
  % fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
toc

%--------------------------------------------------------------------
function [s0,s1,etainf] = svalue(beta)
% Used by the program fsksec
% For a given value of beta, two estimates of
% f''(0), s0 and s1 are computed, together with an
% estimate of etainf. The formulas are given in 
% appendix 3, part 3 of the compendium.
%
p1 = [0.633 -1.68 5.76] ; p2 = [36.76  2.0 5.87];
p3 = [0.125  -0.9 5.463]; betasep = -0.19883768;
% === Compute an estimate of s0 and s1 ===
if(beta <= 1.0)
   s0 = (1.27*(0.2 + beta))^0.56;   % betasep <= beta <= 1.0
   s1 = (1.23*(beta - betasep))^0.54;
else
   s0 = 1.23*(beta^0.454);  % 1.0 < beta <= 1.999
   s1 = (-0.0693*beta + 0.661)*beta + 0.642;
end
% === Compute an estimate of etainf, giving a value
% === of f''(0) around 1.0e-5
if ( beta <= 1.0)
   if ( beta >= 0.0 ) % 0 <= beta <= 1.0
      p = polyval(p1,beta);
      etainf = round(10.0*p)*0.1;
   else
      p = polyval(p2,beta); % betasep <= beta < 0
      etainf = round(10.0*p)*0.1;
   end
else
   p = polyval(p3,beta); % 1 < beta <= 1.999
   etainf = round(10.0*p)*0.1;
end
%----------------------------------------------------------------
function yd = fcnfsk(x,y)
% Used by fsksec
global beta;
yd = zeros(size(y));
yd(1) = y(2);
yd(2) = y(3);
yd(3) = -(y(1)*y(3) + beta*(1 - y(2)^2));
