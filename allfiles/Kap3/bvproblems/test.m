function test
%  Solves the heat-exchanger problem in section 3.3
%  This version uses nested functions

ay = 0.9 ; ai = 0.2;
solinit = bvpinit(linspace(0,1.0,10),@vvinit);

% U(1) = u, U(2) = v
sol = bvp4c(@vvkode,@vvbc,solinit);
% u = sol.y(1,:), v = sol.y(2,:)

FS = 'FontSize';
clf reset

Ty = 100 - 70*sol.y(1,:);
Ti = 100 - 70*sol.y(2,:);
plot(sol.x,Ty,'k',sol.x,Ti,'k-.');
axis([0 1 0 100]);
title('Varmeveksler',FS,14);
xlabel('x',FS,14);
ylabel('T{_i} , T{_y}',FS,14)
legend('T{_y}','T{_i}');
grid
shg

% --------------------------------------------------------------------------
function dUdx = vvkode(x,U)
% the equation for the heat exchanger
% U(1) = u , U(2) = v
dUdx = [ -ay*(U(1) -U(2));  -ai*(U(1) - U(2))];
end

% -------------------------------------------------------------------------
function res = vvbc(U0,U1)
% Boundary conditions for the Falkner-Skan equation
res = [U0(1); U1(2)- 1;];
end
% -------------------------------------------------------------------------
function v = vvinit(x)   
% %  Guessing u = 4*x/7 , v = (x + 6)/7
v = zeros(2,1);
v(1) = 4*x/7;
v(2) = (x + 6)/7;
end
end

% --------------------------------------------------------------------------

