%======= nspring =========
% Used by Nonlinvib
%
function varargout = nspring(t,y,flag)
switch flag
case ''         % Return dy/dt = f(t,y)
   varargout{1} = f(t,y);
case 'events'
   [varargout{1:3}] = events(t,y); % Return [value,isterminal,
                                    %         direction]
otherwise
   error(['Unknown flag ''' flag '''.']);
end
% ------------ Flag is empty -----------

function dydt = f(t,y)
dydt = zeros(size(y));
fac = sqrt(1 + y(1)^2);
dydt(1) = y(2);
dydt(2) = - 2*y(1)*(fac - 1)/fac;

%--------------Flag = 'events' ----------

function [value, isterminal,direction] = events(t,y)
value = y(1);
isterminal = 1;
direction = 0;