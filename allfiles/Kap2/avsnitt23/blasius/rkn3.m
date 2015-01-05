function[t,y] = rkn3(odefun,tspan,y0,h,method,varargin)
% RKN2: Solves a system of ODEs using Euler, three versions of RK2 
% and two versions of RK4. LRH
%
% [t,y] = rkn2(odefun,tspan,y0,h,method) 
% integrates the system of differential equations y' = f(t,y)
% from t0 to tfinal with initial conditions y0.
% The system is coded in the function odefun.
% Function odefun(t,y) must return a column vector
% corresponding to f(t,y).
%
% [t,y] = rkn2(odefun,tspan,y0,h,method,p1,p2,..)
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
%    where y1 = y,y2 = y', y3 = y'' etc.
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
tl = abs(tf - t0); h = abs(h);
n = fix(tl/h); % No. of whole intervals
hlast = mod(tl,h); % Length of last step
nlast = 0;
if hlast > tl*eps*1.0e3
    nlast = 1;
end
% === Test for direction of integration ==
if t0 > tf
    h = -h;
    hlast = -hlast;
end

% === Allocate space ===
y = zeros(n + 1 + nlast,length(y0));
t = zeros(n + 1 + nlast,1);

% === Initial values ===
yvec = y0; y(1,:) = y0';
tt = t0; t(1) = t0; 

switch method
    case {'rk1'}
        % === Forward Euler-scheme ===
        for k = 1: n
            yvec = yvec + h*feval(odefun,tt,yvec);
            y(k+1,:) = yvec';
            tt = t0 + h*k;
            t(k+1) = tt;
        end
        % === Test for a final step ==
        if nlast > 0
            yvec = yvec + hlast*feval(odefun,tt,yvec);
            y(n+2,:) = yvec';
            t(n+2) = tt + hlast;
        end
%         
    case{'rk2GRalston','rk2GMidpoint','rk2'}
        
        if strcmp(method,'rk2')              %% Heun's method
            a1 = 1/2; a2 = 1/2;
            p1 = 1; p2 = 1;
            
        elseif strcmp(method,'rk2GMidpoint') %% Midpoint method
            a1 = 0; a2 = 1;
            p1 = 1/2; p2 = 1/2;
      
        else                               %% Ralston's method
            a1 = 1/3; a2 = 2/3;
            p1 = 3/4; p2 = 3/4;
        
        end
        
        for k = 1: n
            
            k1 = feval(odefun,tt,yvec);
            k2 = feval(odefun,tt+p1*h,yvec+p1*k1*h);
            
            yvec = yvec + (a1*k1 + a2*k2)*h;
            
            tt = t0 + h*k;
            y(k+1,:) = yvec';
            t(k+1) = tt;
        end
        % === Test for a final step ==
        if nlast > 0
            h=hlast;
            k1 = feval(odefun,tt,yvec);
            k2 = feval(odefun,tt+p1*h,yvec+p1*k1*h);
            
            yvec = yvec + (a1*k1 + a2*k2)*h;
            
            tt = tt+h;
            y(n+2,:) = yvec';
            t(n+2) = tt;
        end
        
    case {'rk4'}
        % === Compact RK4-scheme ===
        h2 = 0.5*h;
        for k = 1: n
            y2 = yvec + h2*feval(odefun,tt,yvec);    % 1.stage
            y3 = yvec + h2*feval(odefun,tt + h2,y2); % 2.stage
            y2 = y2 + 2.0*y3;
            y3 = yvec + h*feval(odefun,tt + h2,y3);  % 3.stage
            y2 = y2 + y3;
            yvec = (y2 - yvec + h2*feval(odefun,tt + h,y3))/3; % 4.stage
            tt = t0 + h*k;
            y(k+1,:) = yvec';
            t(k+1) = tt;
        end
        % === Test for a final step ==
        if nlast > 0
            h2 = 0.5*hlast;
            y2 = yvec + h2*feval(odefun,tt,yvec);    % 1.stage
            y3 = yvec + h2*feval(odefun,tt + h2,y2); % 2.stage
            y2 = y2 + 2.0*y3;
            y3 = yvec + hlast*feval(odefun,tt + h2,y3); % 3.stage
            y2 = y2 + y3;
            yvec = (y2 - yvec + h2*feval(odefun,tt + hlast,y3))/3; % 4.stage
            tt = tt + hlast;
            y(n+2,:) = yvec';
            t(n+2) = tt;
        end
    case{'kutta4'}
        for k= 1:n
            k1 = feval(odefun,tt,yvec); 
            k2 = feval(odefun,tt+h/3,yvec+k1*h/3); 
            k3 = feval(odefun,tt+2*h/3,yvec-k1*h/3+h*k2);
            k4 = feval(odefun,tt+h,yvec+k1*h-h*k2+h*k3);
            yvec = yvec + (k1 + 3*k2 + 3*k3 + k4)*h/8;
            tt = t0 + h*k;
            y(k+1,:) = yvec';
            t(k+1) = tt;
        end
        % === Test for a final step ==
        if nlast > 0
            h=hlast;
            k1 = feval(odefun,tt,yvec); 
            k2 = feval(odefun,tt+h/3,yvec+k1*h/3); 
            k3 = feval(odefun,tt+2*h/3,yvec-k1*h/3+h*k2);
            k4 = feval(odefun,tt+h,yvec+k1*h-h*k2+h*k3);
            yvec = yvec + (k1 + 3*k2 + 3*k3 + k4)*h/8;
            tt = tt+h;
            y(n+2,:) = yvec';
            t(n+2) = tt;
        end
        
    otherwise
        error('Unknown value of argument "method" !')
end