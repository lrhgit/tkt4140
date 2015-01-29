% =================== Program cantilev2 ========================
% This is a special version of cantilev1. This version
% focus on the vertical and horisontal deflection of the beam-end.
% (See cantilev1 for a detailed description)
% We have included a (coarse) curve-fitted relation
% between S = T'(0) and the loadfactor alpha^2 = alpha2
%
% Calling function deflect
%
clear all; close all; clc;
FS = 20; LW=3; set(0,'DefaultLineLineWidth',LW,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
global alpha2 

alphstart= 0.0; alphstep = 0.25; alphend = 10.0;
alpha = (alphstart: alphstep: alphend)';
n = length(alpha);
Lh = zeros(n,1); yv = Lh; Lh(1) = 1; % Initialize
S0vec = Lh; S1vec = Lh;
for k = 1:n
   if (alpha(k) > 1)    
      S0vec(k) = 1.15*alpha(k)^0.562;
      S1vec(k) = 0.98*S0vec(k);  
   else 
      S0vec(k) = alpha(k);
      S1vec(k) = 0.95*S0vec(k);
   end   
end  

zspan = [0.0 1.0 ];

% ==== Main Loop ========
%
for k = 2:n  
   % === Compute phi0
   S0 = S0vec(k); S1 = S1vec(k); alpha2 = alpha(k); 
   y0 = [0.0 S0  0.0];
   [z,y] = ode45(@deflect,zspan,y0);
   phi0 = y(end,2);
   % === Initial values to start the iteration  
   itmax = 7; epsi = 1.0e-5; it = 0; dS = 1;
   options = odeset('RelTol',1.0e-5);
   
   % === Start of iteration
   while(abs(dS) > epsi) & (it < itmax)
      it = it + 1;
      y0 = [0.0 S1 0.0];
      [z,y] = ode45(@deflect,zspan,y0,options);
      phi1 = y(end,2);
      dS = -phi1*(S1 - S0)/(phi1 - phi0);
      
      S0 = S1;
      S1 = S1 + dS;
      phi0 = phi1;
      
   end
   yv(k) = y(end,3); 
   % === Compute non-dimensional horisontale length
   Lh(k) = S1/alpha(k);
end
  
% === Plotting yv and Lh as a function of alpha
h=plot(yv,alpha,Lh,alpha,'-.');

%% Improve plot: labels, legends, and title
grid on

xlabel('y');
ylabel('\alpha^2');
h3=legend('\delta','lh');                                 
set(h3,'box','off');
title('Large deflection of cantilever');
