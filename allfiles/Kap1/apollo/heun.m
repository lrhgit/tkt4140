function ynew = heun(fname,x,y,h)
% The 2. order Heun scheme.
% Heun performs one integration step with steplength h.
% Input:
%        x -- The independent variable. Must be
%             given an initial value.
%        y -- Vector containing the current values
%             of the dependent variables
%        h -- The step size
% Output:
%        ynew -- Vector containing the y-values
%                at the end of the integration step.
ynew = zeros(size(y)); 
k1 = feval(fname,x,y);
yp = y + h*k1;
xp = x + h;
k2 = feval(fname,xp,yp);
ynew = y + 0.5*h*(k1 + k2); 

