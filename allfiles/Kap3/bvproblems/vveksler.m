function vveksler
%  Solves the heat exchanger problem in section 3.3

ay = 0.9 ; ai = 0.2;
solinit = bvpinit(linspace(0,1.0,10),@vvinit);
%set 'stats' = 'on' for statistics
options = bvpset('stats','off'); 
sol = bvp4c(@vvkode,@vvbc,solinit,options);
% u = sol.y(1,:), v = sol.y(2,:)

Ty = 100 - 70*sol.y(1,:);
Ti = 100 - 70*sol.y(2,:);
clf
plot(sol.x,Ty,'k',sol.x,Ti,'k-.');
axis([0 1 0 100]);
title('Heat exchanger','FontSize',14);
xlabel('x','FontSize',14);
ylabel('T{_i} , T{_y}','FontSize',14)
legend('T{_y}','T{_i}');
grid
shg
% ----------------------------------------------------
function dUdx = vvkode(x,U)
% the equation for the heat exchanger
% U(1) = u , U(2) = v
dUdx = [ -ay*(U(1) -U(2));  -ai*(U(1) - U(2))];
end
% ----------------------------------------------------
function res = vvbc(U0,U1)
% Boundary conditions : u(0) = 0 , v(1) = 1
res = [U0(1); U1(2)- 1;];
end
% ----------------------------------------------------
function v = vvinit(x)   
% %  Rough guessing
v = zeros(2,1);
v(1) = 0.5;
v(2) = 0.5;
end
end
% ----------------------------------------------------

