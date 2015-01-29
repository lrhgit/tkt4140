% program sor3cont
% Solves example in fig. 7.6 in section 7.3.
% using SOR-iteration.
% This version uses false points 
% along the lines x = 0 and y = 0
% Draw a contor-plot
% h = 0.05;
clear
h = 0.01;
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
epsr = 1.0e-5; % relative iteration error
relres = 1.0; it = 0;
% --- Start iteration ---
while relres > epsr
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
x = [0 : h: 1];
v = [0.1 0.2 0.3 0.4 0.47 0.5 0.53 0.6 0.7 0.8 0.9];
contour(x,x,T',v);
svart = zeros(64,3);
colormap(svart);
%[C,h]=contour(T',v);
%clabel(C,h);
set(gca,'xtick',0:0.2:1.0)
set(gca,'ytick',0:0.2:1.0)
% mesh(x,x,T');
% svart = zeros(64,3);
% colormap(svart);
% set(gca,'ytick',0: 0.2 :1)
% xlabel('x','FontSize',14,'FontWeight','Bold')
% ylabel('y','FontSize',14,'FontWeight','Bold')
% % surf(x,x,T');
% % shading interp
shg

