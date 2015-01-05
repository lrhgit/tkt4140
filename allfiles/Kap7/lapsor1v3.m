% program lapsor1v3
% Solves example in fig. 7.3, section 7.2.
% using SOR-iteration.
% We use (3.13)a as a stopping criterion.
%
% Computes no. of iterations as function of
% omega and the starting-values.
%
% in this version we may select different nets
% by specifing the parameter net.
% net = n + 1
% h(n) = h/2^n , n = 0, 1, 2,..
% starting with h = 0.25 for n = 0 
clear
net = 1;
n = net - 1;
fac = 2^(n+1);
imax = 2*fac + 1; % points in x-direction
jmax = 3*fac + 1; % points in y-direction
startval = 0.0;
fprintf('net = %4.0f \n',net);
fprintf('startvalue = %6.3f \n',startval);
reltol = 1.0e-5; % relative iteration error
for omega = 1.0 :0.01 :1.96
    T = startval*ones(imax,jmax); % initial values
    T(1:imax,jmax) = 100; % boundary values for y = 1.5
    for k = 1:imax % % boundary values for y = 0 
        T(k,1) = 0.0;
    end
    for k = 2:jmax - 1 % % boundary values for x = 0 and 1;
        T(1,k) = 0.0;
        T(imax,k) = 0.0;
    end
    relres = 1.0; it = 0;
    % --- Start iteration ---
    while relres > reltol
        it = it + 1;
        Tsum = 0.0; dTsum = 0.0;
        for i = 2 : imax-1
            for j = 2 : jmax-1
                resid = T(i-1,j) + T(i,j-1) + T(i+1,j) + T(i,j+1)-4*T(i,j);
                dT = 0.25*omega*resid;
                dTsum = dTsum + abs(dT);
                T(i,j) = T(i,j) + dT;
                Tsum = Tsum + abs(T(i,j));
            end
        end
        relres = dTsum/Tsum;
    end
    fprintf('omega = %5.2f  it = %5.0f \n',omega,it);
end

