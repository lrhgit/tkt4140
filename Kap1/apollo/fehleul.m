% Program fehleul
% Sammeligner ode45 med Euler
% Bruker data fra Fehlberg]
clear
global mu lam;
mu = 1/82.45; % 
lam = 1 - mu;
tend = 6.192169331319632; 

% time = [0 tend];
% ystart = [1.2 ; 0; 0; -1.049357509830343]; 
% relfeil = 1.0e-5; absfeil = 1.0e-4*relfeil;
% options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',4);
% [t,y] = ode45(@fcnsonde,time,ystart,options);
% clf;
% plot(y(:,1),y(:,2),'k');

% === Euler og RK4 metode ==
y = [1.2 ; 0; 0; -1.049357509830343]; 
nsteps = 900000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tend/nsteps;
for k = 1:nsteps - 1
    t = k*dt;
    %ynew = euler('fcnsonde',t,y,dt);
    ynew = euler('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
% hold on
 plot(xc,yc,'k'); 
  %hold off
% grid
% axis equal
% xlabel('x','FontSize',14,'FontWeight','Bold')
% ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
