% program lsor4
% Solves example in section 7.3 (7.4)
% using SOR-iteration.
% 
% h = 0.05;
clear
% --- Coarse mesh ---
h = 0.25;
nx = 1/h; ny = nx;
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
Tm = 0.0*ones(imax,jmax); % temperatures
% --- Compute optimal omega ---
ro = cos(pi/nx);
omega = 2/(1 + sqrt(1 - ro^2));
% --- Boundary values ---
for i = 1: imax % along y = 1
    Tm(i,jmax) = 1;
end
for j = 1: jmax - 1 % along x = 1
    Tm(imax,j) = 0;
end
epsr = 1.0e-2; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > epsr
    it = it + 1;
    Tsum = 0.0; dTsum = 0.0;
    for i = 2 : imax-1
        for j = 2: jmax-1
            resid = Tm(i-1,j) + Tm(i,j-1) + Tm(i+1,j) + Tm(i,j+1)-4*Tm(i,j);
            dT = 0.25*omega*resid;
            dTsum = dTsum + abs(dT);
            Tm(i,j) = Tm(i,j) + dT;
            Tsum = Tsum + abs(Tm(i,j));
        end
    end
    Tm(1,1) = 0.5*(Tm(2,1) + Tm(1,2));
    % --- Update boundary values along x = 0 ---
    for i = 2: imax - 1    
        Tm(i,1) = (4*Tm(i,2) - Tm(i,3))/3;
    end
    % --- Update boundary values along y = 0 ---
    for j = 2: jmax - 1    
        Tm(1,j) = (4*Tm(2,j) - Tm(3,j))/3;
    end
    relres = dTsum/Tsum; 
end
Tm
omega
it

% --- Fine mesh ---
h = 0.05;
nx = 1/h; ny = nx;
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
T = 0.5*ones(imax,jmax); % temperatures
% --- Compute optimal omega ---
ro = cos(pi/nx);
omega = 2/(1 + sqrt(1 - ro^2));
% --- Boundary values ---
for i = 1: imax % along y = 1
    T(i,jmax) = 1;
end
for j = 1: jmax - 1 % along x = 1
    T(imax,j) = 0;
end
% --- Interpolation from coarse to fine mesh ---
for i = 1:5
    k = 5*i - 4;
    for j = 1:5
        l = 5*j - 4;
        T(k,l) = Tm(i,j);
    end
end
epsr = 1.0e-4; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > epsr
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
    % --- Update boundary values along x = 0 ---
    for i = 2: imax - 1    
        T(i,1) = (4*T(i,2) - T(i,3))/3;
    end
    % --- Update boundary values along y = 0 ---
    for j = 2: jmax - 1    
        T(1,j) = (4*T(2,j) - T(3,j))/3;
    end
    relres = dTsum/Tsum; 
end
omega
it

