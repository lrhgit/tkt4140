% ================== Program cantilev1 ====================
% Solves for large deflection of a cantilever beam  
% using a shooting technique. 
% The equations are:
%    T''(s) + alpha*cos(T(s))= 0 
% Boundary conditions : T(0) = 0, T'(1) = 0
% Here T denotes the slope angle theta and s the 
% non-dimensional arc length.
% The factor alpha is given by:
%    alpha = P*L^2/EI where P is the vertical
% load at the end, L the length and
% EI the stiffness of the beam.
% The vertical deflection yv is computed from:
%    yv'(s) = sin(T(s)) with yv(0) = 0.
% In the program we put y(1)= T, y(2) = T'
% and y(3) = yv and we use z for the arc length.
% We have to guess S = T'(0) = y(2)(0) in order to satisfy
% phi(S) = 0 with phi(S) = T'(1;S) in this case.
%
% Calling the function deflect
%
clear all; close all; clc;
FS = 20; LW=3; set(0,'DefaultLineLineWidth',LW,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);

global alpha2 
%alpha2  = 3.1393;
alpha2  = 5;
fprintf('\n        Loadfactor alpha^2 = %5.2f \n\n',alpha2  );
zspan = [0.0 1.0 ];
% === From the program cantilevplot we have found
%     suitable initial values S0 og S1 for T'(0)
S0 = 2.0;
S1 = 2.5;

% === Compute phi0
y0 = [0.0 S0  0.0];
[z,y] = ode45(@deflect,zspan,y0);
phi0 = y(end,2);
% === Initial values to start the iteration
itmax = 7; epsi = 1.0e-5; it = 0; dS = 1;
options = odeset('RelTol',1.0e-5);
% === Heading of table
fprintf('        itr.      S          dS\n\n');

% === Start of iteration
while(abs(dS) > epsi) & (it < itmax)
   it = it + 1;
   y0 = [0.0 S1 0.0];
      [z,y] = ode45(@deflect,zspan,y0,options);
   phi1 = y(end,2);
   dS = -phi1*(S1 - S0)/(phi1 - phi0);
   S = S1 + dS;
   S0 = S1;
   S1 = S;
   phi0 = phi1;
   fprintf('%10d %12.6f %12.3e\n',it,S,dS);
end

% === Compute non-dimensional horisontal length
Lh = S/alpha2;
fprintf('\n        Horisontal length Lh = %10.4e \n\n',Lh);
% === Compute a table for theta, theta' and yv 
dz= 0.1; % Print-spacing along the beam 
zspan = (0.0:dz:1.0);
y0 = [0.0 S 0.0];
[z,y] = ode45(@deflect,zspan,y0,options);
fprintf('\n          s       theta     theta''      yvert\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[z y]');

%% plot with labels, legends, and title

h=plot(z,y(:,1));
grid on

xlabel('s');
ylabel('\theta');
%set(hh(3),'box','off');
%title('Shooting with RK4C','FontSize',FS)
