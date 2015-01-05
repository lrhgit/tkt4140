function[t,y] = feuler(odefun,tspan,y0,h)
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
