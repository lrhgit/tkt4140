% program npoisor
% Solves a nonlinear Poisson equation
% on the domain in fig. 7.15, section 7.3
% using SOR-iteration.
% In this version we may select different nets
% by specifying the parameter net.
%
% net = 1 -> h = 0.25, net = 2 -> h = 0.25/2
% giving hn = h/2^(net -1)
clear
net = 2;
h = 0.25;
hn = h/2^(net -1);
nx = 1/hn; ny = nx;
hn2 = hn*hn;
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
% --- Initial values including the boundaries ---
u = zeros(imax,jmax);
%
% --- Compute optimal omega ---
ro = cos(pi/nx);
omega = 2/(1 + sqrt(1 - ro^2));

reltol = 1.0e-5; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > reltol
    it = it + 1;
    usum = 0.0; dusum = 0.0;
    for i = 2 : imax-1
        for j = 2: jmax-1
            fac1 = 4 - 2*hn2*u(i,j) ;
            fac2 = hn2*(u(i,j)^2 + 1);
            resid = u(i-1,j)+u(i,j-1)+u(i+1,j)+u(i,j+1)-4*u(i,j)+fac2;
            du = omega*resid/fac1;
            dusum = dusum + abs(du);
            u(i,j) = u(i,j) + du;
            usum = usum + abs(u(i,j));
        end
    end
    relres = dusum/usum; 
end
u
