% program lap1
% Solves example in fig. 7.3, section 7.2
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

nx = 3;  % # interne nodar i x-retning
ny = 5; % # interne nodar i y-retning

n=nx*ny;  % total number of unknowns
 
xmax=1;  
xmin=0;

ymax=1.5;
ymin=0;

T0=100;

x=linspace(xmin,xmax,nx+2); % nx+2 to account for bondaries
y=linspace(ymin,ymax,ny+2);

d = ones(n,1); % diagonal
b = zeros(n,1); % right hand side

%Impose BCs
b(ny:ny:end)=-T0;

% --- generate A-matrix ---
A = spdiags([d d -4*d d d],[-ny -1 0 1 ny], n,n);
% --- Update A ---

A(ny:ny:end,ny+1:ny:end)=0;
A(ny+1:ny:end,ny:ny:end)=0;
% --- solve system ---


Tvec = A\b;
Temp=reshape(Tvec,ny,nx);
T=zeros(ny+2,nx+2);
T(2:ny+1,2:nx+1) = Temp;
T(end,:) = T0;

figure(1)
%mesh(x(2:nx+1),y(2:ny+1), Temp)
mesh(x,y,T)

xlabel('x')
ylabel('y')


A2 = spdiags([d d -4*d d d],[-nx -1 0 1 nx], n,n);
A2(nx:nx:end,nx+1:nx:end)=0;
A2(nx+1:nx:end,nx:nx:end)=0;

b2 = zeros(n,1); % right hand side
b2(end-nx+1:end) = -T0;


figure(2)
Tvec2=A2\b2;

T2=reshape(Tvec2,nx,ny);
T3=zeros(nx+2,ny+2);
T3(2:nx+1,2:ny+1) = T2;
T3(:,end) = T0;

%mesh(x(2:nx+1),y(2:ny+1), T2')
mesh(x,y,T3')

xlabel('x')
ylabel('y')

%mesh(y(2:ny+1), x(2:nx+1),T2)
