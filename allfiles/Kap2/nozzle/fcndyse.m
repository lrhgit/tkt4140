function fd = fcndyse(x,f)
%--- Plan dyse
% Kalles fra program dyse og plotdyse1
fd = zeros(size(f));
fd(1) = f(2);
fd(2) = f(3);
fd(3) = - (f(1)*f(3) + f(2)^2);