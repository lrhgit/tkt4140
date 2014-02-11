% Program hairer
% Bruker data fra Hairer et al , ref.[15]
clear
global mu lam;
mu = 1/82.45; % Neglisjert månens masse sammelignet med jordas
lam = 1 - mu;
% tendh = 17.065216560157963; % Hairer
tendf = 6.192169331319632; % Fehlberg

% time = [0 tend];
% ystart = [0.994 ; 0; 0;-2.0015851063790825];
% relfeil = 1.0e-5; absfeil = 1.0e-4*relfeil;
% options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',1);
% [t,y] = ode45(@fcnsonde,time,ystart,options);
% clf;
% plot(y(:,1),y(:,2),'k-o');
% === Euler og RK4 metode ==
% y = [0.994 ; 0; 0;-2.0015851063790825]; % Hairer
y = [1.2 ; 0; 0; -1.049357509830343]; % Fehlberg
nsteps = 9000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tendf/nsteps;
for k = 1:nsteps - 1
    t = k*dt;
    %ynew = euler('fcnsonde',t,y,dt);
    ynew = RK4C('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
 plot(xc,yc,'k'); 
% grid
% axis equal
% xlabel('x','FontSize',14,'FontWeight','Bold')
% ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
