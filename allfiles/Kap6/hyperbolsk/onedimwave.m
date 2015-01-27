function [x, t , U] = onedimwave(phi,v,L,A,B,T,N,M,c)
%
% Solves the one-dimensional wave problen
%  d^2*u(x,t)/d^2(t) = c^2*d^2*u(x,t)/d^2(x)
%
%          --- Input ---
% phi(x) : initial wave profile function
% v(x)   : initila wave velocity function
% L      : Length of string
% A(t)   : Height function of left end of string u(0,t) 
% B(t)   : Height function of right end of string u(L,t)
% T      : Final time
% N      : Number of x-grid values
% M      : M = Number of t-grid values
% c      : Speed of wave
%
%          --- Output ---
%  t     : Time grid row vector. Starts at t = 0 and end
%          at t = T with M+2 equally spaced values
%  x     : Space grid row vector
%  U     : (M+2)*(N+2) matrix of solution approximations at 
%          corresponding grid points. Y-grid will correspond to first
%          row indices of U and x-grid values to second column indices
%          of U.
%
h = L/(n+1); % x-step
k = T/(M+1); % t-step
% The Courant C = c*t/h. Stability for  C < 1
U = zeroes(M+2,N+2);
x = 0:h:L; % x-grid
t = 0:k:T; % t-grid
% Assign left and right Dirichlet boundary values.
U(:,1) = feval((A,t))';
U(:,N+2) = feval((B,t)';
% Assign initial time t = 0valuesand next step t = k values
for i = 2:(N+1)
    U(1,i) = feval(phi,x(i));
    mu(i) = 