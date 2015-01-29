% Program fehlberg
% Det avgrensa trelegeme problemet.
% Sammeligner ode45 med RK4 og Eulers metode
% Bruker data fra Fehlberg.
clear
global mu lam;
mu = 1/82.45; % 
lam = 1 - mu;
tend = 6.192169331319632; 

% === Eulers metode ==
y = [1.2 ; 0; 0; -1.049357509830343]; 
nsteps = 900000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tend/nsteps;
tic
for k = 1:nsteps - 1
    t = k*dt;
    ynew = euler('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
toc
clf;
plot(xc,yc,'k--'); 

% === RK4C  ==
y = [1.2 ; 0; 0; -1.049357509830343]; 
nsteps = 3000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tend/nsteps;
tic
for k = 1:nsteps - 1
    t = k*dt;
    ynew = RK4C('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
toc
hold on
plot(xc,yc,'k-.'); 

% === ode45  ==
time = [0 tend];
ystart = [1.2 ; 0; 0; -1.049357509830343]; 
relfeil = 1.0e-4; absfeil = 1.0e-4*relfeil;
options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',4);
tic
[t,y] = ode45(@fcnsonde,time,ystart,options);
toc
plot(y(:,1),y(:,2),'k');
hold off
grid
axis equal
xlabel('x','FontSize',14,'FontWeight','Bold')
ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
