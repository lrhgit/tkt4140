function ynew = rk4c(fname,x,y,h)
% The classical RK4-scheme.
% RK4C performs one integration step with steplength h.
% Input:
%        x -- The independent variable. Must be
%             given an initial value.
%        y -- Vector containing the current values
%             of the dependent variables
%        h -- The step size
% Output:
%        ynew -- Vector containing the y-values
%                at the end of the integration step.

k1 = zeros(size(y)); k2 = k1; k3 = k1;
k4 = k1; ytemp = k1; ynew = k1; 
hhalf = 0.5*h;

k1 = h*feval(fname,x,y);           % 1. stage
ytemp = y + k1/2;
k2 = h*feval(fname,x+hhalf,ytemp); % 2. stage
ytemp = y + k2/2;
k3 = h*feval(fname,x+hhalf,ytemp); % 3. stage
ytemp = y + k3;
k4 = h*feval(fname,x+h,ytemp);     % 4. stage
ynew = y + (k1 + 2*k2 + 2*k3 + k4 )/6.0;
