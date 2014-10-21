function[t,y] = heun(odefun,tspan,y0,h)
t0 = tspan(1); tf = tspan(2);
tl = abs(tf - t0); h = abs(h);
n = fix(tl/h); % No. of whole intervals
hlast = mod(tl,h); % Length of last step
nlast = 0;
if hlast > tl*eps*100
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

% === Heun-scheme ===
for k = 1: n
    val = feval(odefun,tt,yvec);
    ypred = yvec + h*val;
    tt = t0 + h*k;
    yvec = yvec + 0.5*h*(val + feval(odefun,tt,ypred));
    y(k+1,:) = yvec';
    t(k+1) = tt;
end

% === Test for a final step ==
if nlast > 0
    val = feval(odefun,tt,yvec);
    ypred = yvec + hlast*val;
    tt = tt + hlast;
    yvec = yvec + 0.5*hlast*(val + feval(odefun,tt,ypred));
    y(n+2,:) = yvec';
    t(n+2) = tt;
end
