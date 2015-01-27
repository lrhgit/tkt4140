function[t,y] = rkn(odefun,tspan,y0,h,method,varargin)
% RKN: Solves a system of ODEs using Euler, Heun and RK4.
%
% [t,y] = rkn(odefun,tspan,y0,h,method) 
% integrates the system of differential equations y' = f(t,y)
% from t0 to tfinal with initial conditions y0.
% The system is coded in the function odefun.
% Function odefun(t,y) must return a column vector
% corresponding to f(t,y).
%
% [t,y] = rkn(odefun,tspan,y0,h,method,p1,p2,..)
% passes the additional parameters p1,p2,.. to the function odefun
% as odefun(t,y,p1,p2,...)
% Note : The parameters p1,p2,.. must follow after the argument "method".
%
% === Input === 
%
% tspan = [t0 tfinal]: If t0 > tfinal, the integration goes
%                       from right to left.
% y0 is a column vector containing the initial values
%    for the dependent variables y1,y2,i.e : y1(t0),y2(t0),..
% h is the steplength: Any positive value < abs(t0 -tfinal)
% method : A string with the following values:
%          'rk1' : Forward Euler
%          'rk2' : Heun 
%          'rk4' : RK4
%
% === Output ===
% 
% t: a column vector containing the values of the
%    independent variable.
% y: a matrix where the first column contains the
%      values of y1, the next column the values of y2 etc.
%
t0 = tspan(1); tf = tspan(2);
h = abs(h); ht = h;
if t0 > tf
    ht = - ht;
end
t = (t0 : ht: tf)';
tl = abs(tf - t0);
if mod(tl,h) > tl*1.0e-8
    t = [t; tf];
end
n = length(t) - 1;
hvec = diff(t);

% === Allocate space ===
y = zeros(n + 1 ,length(y0));

% === Initial values ===
yvec = y0; y(1,:) = y0'; 

switch method
    case {'rk1'}
        % === Forward Euler-scheme ===
        for k = 1: n
            h = hvec(k);
            tt = t(k);
            yvec = yvec + h*feval(odefun,tt,yvec,varargin{:});
            y(k+1,:) = yvec';
        end
    case {'rk2'}
        % === Heun-scheme ===
        for k = 1: n
            h = hvec(k);
            tt = t(k);
            val = feval(odefun,tt,yvec,varargin{:});
            ypred = yvec + h*val;
            tt = t(k+1);
            yvec = yvec + 0.5*h*(val + feval(odefun,tt,ypred,varargin{:}));
            y(k+1,:) = yvec';
        end
    case {'rk4'}
        % === Compact RK4-scheme ===
        for k = 1: n
            h = hvec(k);
            tt = t(k);
            h2 = 0.5*h;
            y2 = yvec + h2*feval(odefun,tt,yvec,varargin{:});    % 1.stage
            y3 = yvec + h2*feval(odefun,tt + h2,y2,varargin{:}); % 2.stage
            y2 = y2 + 2.0*y3;
            y3 = yvec + h*feval(odefun,tt + h2,y3,varargin{:});  % 3.stage
            y2 = y2 + y3;
            yvec = (y2 - yvec + h2*feval(odefun,tt + h,y3,varargin{:}))/3; % 4.stage
            y(k+1,:) = yvec';
        end
    otherwise
        error('Unknown value of argument "method" !')
end