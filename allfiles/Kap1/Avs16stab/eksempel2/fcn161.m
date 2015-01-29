function varargout = fcn161(t,y,flag)
switch flag
case ''
   varargout{1} = f(t,y);
case 'events'
   [varargout{1:3}] = events(t,y); 
case 'jacobian'
   varargout{1} = jacobian(t,y);  
otherwise
   error(['Unknown flag ''' flag '''.']);
end
% -----------------------------------------------
function dydt = f(t,y)   
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -129600.0*y(2) + 98696.0*y(1)+ 9869.6;
% -----------------------------------------------
function [value,isterminal,direction] = events(t,y)
% Finn tiden når partikkelen treffer ytre sylinder
value = [0.01 - y(1); y(2)];
isterminal = [1; 0];
direction = [0 ; 0];
% ------------------------------------------------
function dfdy = jacobian(t,y)
dfdy = [0 1 ; 98696.0 -129600.0];
% -----------------------------------------------