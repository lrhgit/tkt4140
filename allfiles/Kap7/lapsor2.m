% program lapsor2
% Solves example in fig. 7.6, section 7.2.2
% using SOR-iteration and with false points 
% along the lines x = 0 and y = 0
% In this version we may select different nets
% by specifying the parameter net.
%
% net = 1 -> h = 0.25, net = 2 -> h = 0.25/2
% giving hn = h/2^(net -1)
clear
net = 1;
h = 0.25;
hn = h/2^(net -1);
nx = 1/hn; ny = nx;
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
T = 0.501*ones(imax,jmax); % temperatures
% --- Compute optimal omega ---
ro = cos(pi/nx);
omega = 2/(1 + sqrt(1 - ro^2));

T(1:imax,jmax) = 1;  % boundary values along y = 1
T(imax,1:jmax-1) = 0;% boundary values along x = 1

reltol = 1.0e-5; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > reltol
    it = it + 1;
    Tsum = 0.0; dTsum = 0.0;
    % --- boundary values along y = 0 ---
    for i = 2: imax - 1 
        resid = 2*T(i,2) + T(i-1,1) + T(i+1,1) - 4*T(i,1);
        dT = 0.25*omega*resid;
        dTsum = dTsum + abs(dT);
        T(i,1) = T(i,1) + dT;
        Tsum = Tsum + abs(T(i,1));
    end
    % --- boundary values along x = 0 ---
    for j = 2: jmax - 1    
        resid = 2*T(2,j) + T(1,j-1) + T(1,j+1) - 4*T(1,j);
        dT = 0.25*omega*resid;
        dTsum = dTsum + abs(dT);
        T(1,j) = T(1,j) + dT;
        Tsum = Tsum + abs(T(1,j));
    end
    for i = 2 : imax-1
        for j = 2: jmax-1
            resid = T(i-1,j) + T(i,j-1) + T(i+1,j) + T(i,j+1)-4*T(i,j);
            dT = 0.25*omega*resid;
            dTsum = dTsum + abs(dT);
            T(i,j) = T(i,j) + dT;
            Tsum = Tsum + abs(T(i,j));
        end
    end
    T(1,1) = 0.5*(T(2,1) + T(1,2));
    relres = dTsum/Tsum; 
end

