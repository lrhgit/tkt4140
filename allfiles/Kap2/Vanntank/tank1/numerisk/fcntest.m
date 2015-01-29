function varargout = fcntank(x,y,flag,beta4)
% Brukt i programmet tank
switch flag
case ''
   varargout{1} = f(x,y,beta4); 
otherwise
   error(['Unknown flag ''' flag '''.']);
end
% ---------------------------------------
function dydx = f(x,y,beta4)   
dydx = zeros(size(y));
dydx(1) = y(2);
dydx(2) = y(3);
dydx(3) = y(4);
dydx(4) = -4*beta4*(y(1) + 1 - x);
% ---------------------------------------
