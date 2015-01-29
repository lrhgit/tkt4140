% program lapsor1v3
% Solves example in fig. 7.3, section 7.2.
% using SOR-iteration.
% We use (3.13)a as a stopping criterion.
% This version uses different values
% for the space-increment h
% h(n) = h/2^n , n = 0, 1, 2,..
% starting with h = 0.25 for n = 0 
clear
n = 3;
fac = 2^(n+1);
imax = 2*fac + 1; % points in x-direction
jmax = 3*fac + 1; % points in y-direction

reltol = 1.0e-5; % relative iteration error
for omega = 1.0 :0.05 :1.95
    T = zeros(imax,jmax); % temperatures
    T(1:imax,jmax) = 100; % boundary values
    relres = 1.0; it = 0;
    % --- Start iteration ---
    while relres > reltol
        it = it + 1;
        Tmax = 0.0; dTmax = 0.0;
        for i = 2 : imax-1
            for j = 2: jmax-1
                resid = T(i-1,j) + T(i,j-1) + T(i+1,j) + T(i,j+1)-4*T(i,j);
                dT = 0.25*omega*resid;
                dTmax = max(dTmax, abs(dT));
                T(i,j) = T(i,j) + dT;
                Tmax = max(Tmax, abs(T(i,j)));
            end
        end
        relres = dTmax/Tmax;
    end
    fprintf('omega = %5.2f  it = %5.0f \n',omega,it);
end

