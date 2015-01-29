function startupv3 
% ====== Startup of a flow in a pipe (Szymanski's problem) ===
% See section 5.6 in the compendium.
% Because the left hand side of the system of equations is independent
% of time, we use the modified version of tdma given in 
% eq. (13), appendix 9. 
% This version is two to three times
% faster than the original startup
%
% The equation of w is given by :
%
%     dw/dt = w'' + w'/r , w = w(r,t), 0 < r < 1
%
% Velocity-field : u(r,t) = us(r) - w(r,t)
%                  where us = 1 - r^2
% Boundary conditions : w(1) = w(-1) = 0, w'(0) = 0 
% Initial condition : w(r,0) = 1- r^2 = us
% Timestep: dt = D*h^2 , h = dr
%
% Using nested functions (Matlab ver. 7.x)
% theta = 0:   FTCS-scheme
% theta = 1/2: Crank-Nicolson 
% theta = 1:   Laasonen
%
% Output of u. No  plotting.
% Local functions : solut  
% External functions : fcnwa (analytical solution)
%============================================================
clear
N = 50; % No of parts
h = 1/N; % Length of r-step
D = 0.4; % Numerical diffusion number
dt = D*h^2;
nmax = 500; % No. of time-steps
tid = nmax*dt;
disp('      *****************************************');
disp('      *    Impulsive start of pipeflow        *');
disp('      *    theta = 0:   FTCS-scheme           *');
disp('      *    theta = 1:   Laasonen              *');
disp('      *    theta = 1/2: Crank-Nicolson        *');
disp('      *****************************************');
disp('');
fprintf('\n No. of time-steps...........  %4.0f\n',nmax);
fprintf(' r-step length.................. %7.3f\n',h);
fprintf(' Numerisk diffusjonstall........ %7.3f\n',D);
fprintf(' Timestep....................... %12.3e\n',dt);
fprintf(' Elapsed time................... %12.3e\n\n',tid);
%
% --- Allocating vectors 
a = zeros(N,1); b = a; c = a; d = a;
wnew = a;
wa = zeros(N + 1,1);
wold = wa;

r = (0:h:1)';

% === theta = 0 : FTCS ===
theta = 0; wold = 1 - r.^2;
w1 = solut(theta);

% === theta = 1: Laasonen ===
theta = 1; wold = 1 - r.^2;
w2 = solut(theta);

% === theta = 1/2 : Crank-Nicolson ===
theta = 0.5; wold = 1 - r.^2;
w3 = solut(theta);

% === Analytical solution ===
for l = 1: length(r)
    x = r(l);
    wa(l) = fcnwa(x,tid);
end
%=== Output of av u1, u2 , u3 og ua ===
w1 = [w1; 0]; w2 = [w2; 0]; w3 = [w3; 0];
w1 = 1 - r.^2 - w1; % u for FTCS
w2 = 1 - r.^2 - w2; % u for Laasonen
w3 = 1 - r.^2 - w3; % u Crank-Nicolson
wa = 1 - r.^2 - wa; % u analytical
fprintf('      r       u(ftcs)      u(laasonen)   u(crank-n)    analyt  \n\n')
fprintf('  %7.3f  %12.4e  %12.4e  %12.4e  %12.4e \n',[r w1 w2 w3 wa]');

% === Local function solut solves the theta-scheme ===
%
function vec = solut(theta)
tmp1 = D*(1 - theta);
for j = 2:N
    fac = 0.5/(j - 1);
    a(j) = -D*theta*(1 - fac);
    b(j) = (1 + 2*D*theta) ;
    c(j) = - D*theta*(1 + fac);
end
b(1) = 1 + 4*D*theta;
c(1) = - 4*D*theta;
%=== Elimination,eq.(13a),appendix 9
c(1) = c(1)/b(1);
for j = 2 : N
    b(j) = b(j)- a(j)*c(j-1);
    c(j) = c(j)/b(j);
end
for n = 1: nmax
    % === Compute the right hand side 
    for j = 2:N       
        fac= 0.5/(j - 1);
        d(j) = tmp1*((1 - fac)*wold(j-1)+ (1 + fac)*wold(j+1))...
               +(1 - 2*tmp1)*wold(j);
    end
    d(1) = wold(1) + 4*tmp1*(wold(2) - wold(1));
    % === Backsubstitution , eq. 13b and c 
    d(1) = d(1)/b(1);
    for j = 2:N
        d(j) = (d(j) - a(j)*d(j-1))/b(j);
    end
    wnew(N) = d(N);
    for j = N-1:-1 :1
        wnew(j) = d(j) - c(j)*wnew(j+1);
    end
    wold = [wnew;0]; % Next time-step
end
vec = wnew;
end % function solut
end % function startupv3
