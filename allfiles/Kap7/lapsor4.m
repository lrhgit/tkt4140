% program lapsor4
% Solves example in fig. 7.6, section 7.3.
% using SOR-iteration.
% This version uses 2. order forward
% differences along the lines x = 0 and y = 0
% h = 0.05;
clear
h = 0.05;
nx = 1/h; ny = nx;
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
T = 1.0*ones(imax,jmax); % temperatures
% --- Compute optimal omega ---
ro = cos(pi/nx);
omega = 2/(1 + sqrt(1 - ro^2));

T(1:imax,jmax) = 1;% boundary values along y = 1
T(imax,1:jmax-1) = 0; % boundary values along x = 1

reltol = 1.0e-4; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > reltol
    it = it + 1;
    Tsum = 0.0; dTsum = 0.0;
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
    % --- Update boundary values along y = 0 ---    
    T(2:imax-1 ,1) = (4*T(2:imax-1,2) - T(2:imax-1,3))/3;
    % --- Update boundary values along x = 0 ---   
    T(1,2:jmax-1 ) = (4*T(2,2:jmax-1) - T(3,2:jmax-1))/3;
    relres = dTsum/Tsum; 
end
it

