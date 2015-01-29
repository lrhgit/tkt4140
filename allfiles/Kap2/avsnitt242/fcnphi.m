 function phi = fncphi(r,s)
 % Solves the Blasius equation for two uniform mixing layers.
 % This is the function for the upper layer with eta in [0 , etamax]
 global  odetol u2u1 etamax ksimax
 f0 = [ 0; r; s];
 etaspan = [0 etamax];
 options = odeset('RelTol',odetol);
 [eta,f] = ode45(@fcn,etaspan,f0,options); 
 phi = f(end,2) - 1 ;
 function dy = fcn(x,y)
 dy = zeros(size(y));
 dy(1) = y(2);
 dy(2) = y(3);
 dy(3) = -y(1)*y(3); 
