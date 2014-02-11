function col143
%  Solves the problem page 143 and 182
% in L.Collatz : "The Numerical Treatment of
% Differential Equations ", 3. utgave, Springer 1960
%
% Differential equation :
% m''(x) + (1 + x^2)*m(x) + 1 = 0 
% Boundary conditions: m(1) = m(-1) = 0
%
% m(x) is the moment-distribution of a column-beam 
% with av parabolic varying bending stiffness.
%
a1 = 0.92925; a2 = - 0.05115;
solinit = bvpinit(linspace(-1.0,1.0,20),@colinit);
%set 'stats' = 'on' for statistics
options = bvpset('stats','off'); 
sol = bvp4c(@colkode,@colbc,solinit,options);
% y = sol.y(1,:), y' = sol.y(2,:)
clf
plot(sol.x,sol.y(1,:)','k');
% axis([0 1 0 100]);
FS = 'FontSize';
title('Moment distribution in a parabolic columnbeam. ',FS,12);
xlabel('x',FS,14);
ylabel('m',FS,14,'Rotation',0)
grid
shg
% ----------------------------------------------------
function dmdx = colkode(x,m)
% m(1) = y , m(2) = y'
dmdx = [m(2);  -((1 + x^2)*m(1) + 1)];
end
% ----------------------------------------------------
function res = colbc(mm1,mp1)
% Boundary conditions : m(-1) = 0 , m(1) = 0
res = [mm1(1); mp1(1)];
end
% ----------------------------------------------------
function m = colinit(x)   
m = zeros(2,1);
m(2) = 2*x*((1 -2*x^2)*a2 - a1);
m(1) = (1 - x^2)*(a1 + a2*x^2);
end
end
% ----------------------------------------------------

