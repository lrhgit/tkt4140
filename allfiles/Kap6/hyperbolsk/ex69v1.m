% program ex69v1
% Computes the numerical solution of the example in
% section 6.9, equation 9.3b
%
% Analytical solution:
% u(x,t) = cos(t)*sin(x)
% der 0 <= x <= pi
%
clear; close
x = linspace(0,pi,101)';
M = length(x) - 1;
C = 1.1; % The Courant-number
C2 = C*C;
dx = pi/M;
dt = C*dx;
ustart = sin(x);% Initial values 
unp1 = zeros(M+1,1);un = unp1; unm1 = unp1;
J = 2:M ; % Index-vector
nmax = 40; % No. of time steps
% First time step, eq.9.5
unm1 = ustart;
un(J) = C2*(ustart(J-1) + ustart(J+1))*0.5 + (1 - C2)*ustart(J);
% Further computation using eq.9.4
for n = 2 : nmax 
    t = n*dt;
    unp1(J) = C2*(un(J-1) + un(J+1)) + 2*(1 - C2)*un(J) - unm1(J);
    unm1 = un;
    un = unp1;
end
fprintf(' t = %6.4f \n\n',t);
ua = cos(t)*sin(x); % Analytical
plot(x,unp1,'k',x,ua,'k')
grid
axis([0 pi 0 0.3])
FS = 'FontSize';
xlabel('x',FS,14)
ylabel('u',FS,14,'Rotation',0)
title('u(x,t) = cos(t)sin(x)',FS,14)

        
