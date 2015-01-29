%%% 
% A program for solution of the van der Pol equation
% Different versions of rk2 may be compared by selecting rkfamily='rk2'
% and two versions of rk4 be compared by selecting rkfamily='rk4'
% 
% For for info on the van der Pol equation:   
%                   http://en.wikipedia.org/wiki/Van_der_Pol_oscillator
%
clear all; close all; clc;
FS = 20; set(0,'DefaultLineLineWidth',3,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
%% Select family of methods to compare: rk2 for RK second order 
%                                       rk4 for RK fourh order
rkfamily='rk2';
rkfamily='rk4';

%% Select timestep
dt = 0.125;

tspan = [0, 40];
y0 = [2; 0];
Mu = 5;

rhsf='vdpoolf';
%% Euler method
[t1,ye]  = rkn2(rhsf,tspan,y0,dt,'rk1',Mu);
hf(1)=figure(1);

switch rkfamily
    
    case{'rk2'} %% RK second order
        [t1,yh]  = rkn2(rhsf,tspan,y0,dt,'rk2',Mu);
        [t1,yra] = rkn2(rhsf,tspan,y0,dt,'rk2GRalston',Mu);
        [t1,yrm] = rkn2(rhsf,tspan,y0,dt,'rk2GMidpoint',Mu);
        h=plot(t1,yh(:,1),t1,yrm(:,1),':',t1,yra(:,1),t1,yra(:,1),':');
       
        hh(3,:)=legend('Heun','Midpoint','Ralston');

    case{'rk4'} %% Higher order methods
        [t1,yrk] = rkn2(rhsf,tspan,y0,dt,'rk4',Mu);
        [t1,yk]  = rkn2(rhsf,tspan,y0,dt,'kutta4',Mu);
        
        ode = @(t,y) vanderpoldemo(t,y,Mu);
        [t,y] = ode45(ode, tspan, y0);

        h=plot(t,y(:,1),t1,yrk(:,1),t1,yk(:,1));
        hh(3,:)=legend('ode45','RK4','Kutta4');
    
    otherwise
        error('Unknown value family of methods!')
end

%% Labels and title
hh(1,:)=xlabel('t');
hh(2,:)=ylabel('solution y');

% Title 
titleString=sprintf('van der Pol Equation, \\mu = %5.2f',Mu);
title(titleString)

%% Plot of the phase portrait 
hf(2)=figure(2);

switch rkfamily
    case{'rk2'}
        h2=plot(yh(:,1),yh(:,2),yrm(:,1),yrm(:,2),yra(:,1),yra(:,2));
        hh(6,:)=legend('Heun','Midpoint','Ralston' );

    case{'rk4'}
        h2=plot(y(:,1),y(:,2),yrk(:,1),yrk(:,2),yk(:,1),yk(:,2));
        hh(6,:)=legend('ode45','RK4','Kutta4');

    otherwise
        error('Unknown value family of methods!')
end

set(h2(:),'linewidth',2);
hh(4,:)=xlabel('position');
hh(5,:)=ylabel('velocity');

set(hh(6),'box','off');
Position=get(hf(2),'Position');
dp=100;
Position(1)=Position(1)+dp;
set(hf(2),'Position',Position);
% Title 
titleString=sprintf('Phase portrait van der Pol Equation, \\mu = %5.2f',Mu);
title(titleString)


