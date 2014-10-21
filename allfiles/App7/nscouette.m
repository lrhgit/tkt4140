% program nscouette
%
% Computes the numerical solution of the non-steady
% Couette flow (or heat conduction problem) in
% section 5.2 of the compendium using the FTCS-scheme.
% Equation : du/dt = u''(y,t), 0 <= y <= 1
% Initial values : u(y,t) = 0 , t < 0
% Boundary values : u(0,t) = 1 , u(1,t) = 0
% The relative error (in %) between
% the numerical and analytical is also computed.
% The analytical solution u(y,t) is computed 
% using the function fcnu.
% r = dt/dy^2 and dy is given
% dt is computed from dt = dy^2*r
% Computation from t = 0 to t = tf
%
clear
tf = 0.03;   % final value of t
r = 0.3;     %  numerical Fourier number
dy = 0.1;    % space-step
dt = dy^2*r; % time step

nf = round(tf/dt); % no. of time steps
yv = (0 : dy : 1)' ;
m = length(yv);
u = zeros(m,1);     % initial values
ua = u; err = u;    % allocate vectors
u(1) = 1; u(m) = 0; % boundary values

for k = 1 : nf
    t = k*dt;
    a = u(1);
    b = u(2);
    for j = 2 : m - 1 
        u(j) = r*(a + u(j+1)) + (1 - 2*r)*b;
        a = b;
        b = u(j+1);
    end   
end 
fprintf('t = %7.4f \n\n',t);
fprintf('%5.2f %10.4f \n',[yv u]');
% --- Relative error in % ---
for j = 2: m -1
    y = yv(j);
    ua(j) = fcnu(y,t);
    err(j) = (u(j) - ua(j))/ua(j)*100;
end
fprintf('\n');
fprintf('%5.2f %10.2f \n',[yv err]');
