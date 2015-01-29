% program difftest
clear;
% fprintf('        Kinematical viscosity . nu   = %10.3e m^2/s \n',nu );
% fprintf('        Density of fluid ...... rof  = %10.3e kg/m^3 \n',rof);
% fprintf('        Density of sphere ..... ros  = %10.3e kg/m^3 \n',ros);
% fprintf('        Diameter of sphere .... d    = %10.3e m \n',d);
% fprintf('        Initial velocity....... v0   = %10.3e m/s \n',v0);
% fprintf('        x-comp. of fluid vel... vfx  = %10.3e m/s \n',vfx);
% fprintf('        y-comp. of fluid vel... vfy  = %10.3e m/s \n',vfy);
% fprintf('        Angle of elevation..... alf  = %10.3e deg. \n',alf);

xint = [0 3]; % x-interval
y0 = [-8.0 ]; % Initial values
options = odeset('RelTol',1.0e-5);
[x,y] = ode45(@fcntest,xint,y0,options);
%disp([te ye ie]);
% === Plotting the trajectory ===
% FS = 'FontSize'; FW = 'FontWeight';
plot(x,y,'k-o');
grid
axis equal
% ylim([0 Inf])
% xlabel('x(m)',FS,14,FW,'Bold')
% ylabel('y(m)',FS,14,FW,'Bold')
%
