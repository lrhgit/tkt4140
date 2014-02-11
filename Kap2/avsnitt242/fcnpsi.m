 function psi = fcnpsi(r,s)
% Solves the Blasius equation for two uniform mixing layers.
% The function is for the lower layer with ksi in [0 , ksimax]
 global odetol u2u1 etamax ksimax
 
 f0 = [0; -r; s];
 ksispan = [0 ksimax];
 options = odeset('RelTol',odetol);
 [ksi,f] = ode45(@fcn,ksispan,f0,options);
 psi = f(end,2) + u2u1;
 function dy = fcn(x,y)
 dy = zeros(size(y));
 dy(1) = y(2);
 dy(2) = y(3);
 dy(3) = y(1)*y(3); 

