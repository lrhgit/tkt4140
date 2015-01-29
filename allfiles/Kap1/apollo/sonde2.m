% Program sonde2
% Tester forskjellige baner
% med Robert Newtons verdier
% Merk at Newtons bane nr.7 ikke er
% med,slik at bane = 7 gir Newtons bane 8 osv.
clear; close
clear global mu lam;
global mu lam;
mu = 1/82.45; % R. Newton og Fehlberg
lam = 1 - mu;
%mu = 1/82.2845; % Nyere data fra NASA

% Data fra artikkel av Robert E. Newton
% fprintf('=== Newton-data === \n\n')
start = zeros(10,5);
start(1:10,1) = [ 1.0 1.05 1.1 1.15 1.2 1.0 1.05 1.1 1.15 1.2]';
start(1:10,4) = -[2.3314 1.7739 1.6604 1.5842 1.498 1.5364 0.8475 0.8303 0.9124 1.049]';
start(1:10,5) =  [7.8925 6.5227 6.3841 6.3335 6.302 5.4292 5.7574 5.9750 6.1081 6.192]';

banenr = 10; % 1 <= bane <= 10

time = [0 start(banenr,5)];
ystart = [start(banenr,1) ;start(banenr,2);start(banenr,3);start(banenr,4)];
%relfeil = 1.0e-6; absfeil = 1.0e-4*relfeil;
relfeil = 1.0e-4; absfeil = 1.0e-2*relfeil;
options = odeset('RelTol',relfeil,'AbsTol',absfeil,'Refine',4);
[t,y] = ode45(@fcnsonde,time,ystart,options);
plot(y(:,1),y(:,2),'k-o');
grid
axis equal
xlabel('x','FontSize',14,'FontWeight','Bold')
ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
