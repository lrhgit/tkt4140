function ynew = RK4(fname,x,y,h)
% Compact version of the classical RK4-scheme.
% RK4 performs one integration step with steplength h.
% Input:
%        x -- The independent variable. Must be
%             given an initial value.
%        y -- Vector containing the current values
%             of the dependent variables
%        h -- The step size
% Output:
%        ynew -- Vector containing the y-values
%                at the end of the integration step.
% 
y2 = zeros(size(y)); y3 = y2;
hhalf = 0.5*h;
y2 = y + hhalf*feval(fname,x,y); % 1. stage
y3 = y + hhalf*feval(fname,x+hhalf,y2);% 2. stage
y2 = y2 + 2.0*y3; 
y3 = y + h*feval(fname,x+hhalf,y3); % 3. stage
y2 = y2 + y3;
ynew = (y2 - y + hhalf*feval(fname,x+h,y3))/3.0; % 4. stage


