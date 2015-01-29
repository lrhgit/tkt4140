function ynew = euler(fname,x,y,h)
% The basic Euler-scheme.
% euler performs one integration step with steplength h.
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
ynew = y + h*feval(fname,x,y); 
