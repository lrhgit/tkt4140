% Program hairer
% Det avgrensa trelegeme problemet.
% Sammeligner ode45 med RK4 og Eulers metode
% Bruker data fra Hairer et al, ref. [15]. side 130
clear
global mu lam;
mu = 1/81.45; % Neglisjert månens masse i forhold til jordas.
lam = 1 - mu;
tend = 17.065216560157963; 
clf;

% % === Heuns metode ==
y = [0.994 ; 0; 0; -2.0015851063790825]; 
nsteps = 50000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tend/nsteps;
for k = 1:nsteps - 1
    t = k*dt;
    ynew = heun('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
plot(xc,yc,'k--'); 
hold on
% === RK4C  ==
y = [0.994 ; 0; 0; -2.0015851063790825]; 
nsteps = 10000 ; % Antall skritt
xc = zeros(nsteps,1); yc = xc; 
xc(1) = y(1); yc(1) = y(2);
dt = tend/nsteps;
for k = 1:nsteps - 1
    t = k*dt;
    ynew = RK4C('fcnsonde',t,y,dt);
    y = ynew;
    xc(k+1) = y(1);
    yc(k+1) = y(2);
end
 plot(xc,yc,'k-.'); 

% === ode45  ==
time = [0 tend];
ystart = [0.994 ; 0; 0; -2.0015851063790825]; 
relfeil = 1.0e-4; absfeil = 1.0e-4*relfeil;
options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',4);
[t,y] = ode45(@fcnsonde,time,ystart,options);
plot(y(:,1),y(:,2),'k');
hold off
grid
axis equal
xlabel('x','FontSize',14,'FontWeight','Bold')
ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
legend('Heun: 50000 skritt','RK4: 10000 skritt ','ODE45 : 92 skritt ',2);
