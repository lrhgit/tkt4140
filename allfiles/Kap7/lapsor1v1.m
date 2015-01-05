% program lapsor1v1
% Solves example in fig. 7.3, section 7.2
% using SOR-iteration
clear
nx = 4 ; % parts in x-direction
ny = 6;  % parts in y-direction
imax = nx + 1; % points in x-direction
jmax = ny + 1; % points in y-direction
T = zeros(imax,jmax); % temperature-matrix

T(1:imax,jmax) = 100; % boundary values 
omega = 1.5; % relaxation factor
% --- Start iteration ---
for it = 1: 20
    for i = 2 : imax-1
        for j = 2: jmax-1
            resid = (T(i-1,j)+T(i,j-1) + T(i+1,j) + T(i,j+1) - 4*T(i,j));
            dT = 0.25*omega*resid;
            T(i,j) = T(i,j) + dT;
        end
    end
end
