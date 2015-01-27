% ===================== Program jhamel ============================
% Solution of Jeffrey-Hamel-flow, appendix 5,
% using a shooting technique with the secant-method as a root-finder.
% The Jeffrey-Hamel equation is given by :
%  f'''(eta) + 2*Re*alpha*f(eta)*f'(eta) + 4*alpha^2*f'(eta) = 0
%  f(0) = 1, f'(0) = 0, f(1)=0
% In the  program we put y(1)= f, y(2) = f' and y(3) = f'' 
% 
% s = f''(0) = y(3)(0)
%
% Using function fcnjh
% =====================================================================
clear
global Re alpha
Rea = 0.5;
alphag = 10 ; % Half of the opening-angle in degrees.
alpha = alphag*pi/180; % Converting to radianes.
Re = Rea/alpha;
etaspan = [0 1]; % The integration-interval
% We have found suitable starting-values s0 and s1
% from the program plotjh.
s0 = -4.5;
s1 = -3.5;
% ---- Compute phi0 ----
f0 = [1.0 0.0  s0];
[eta,f] = ode45(@fcnjh,etaspan,f0);
phi0 = f(end,1);
% Starting values for the iteration
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
% Print heading of table
fprintf('        itr.      s          ds\n\n');
% -------- Start of iteration -------------------
while(abs(ds) > epsi) & (it < itmax)
   it = it + 1;
   f0 = [1.0 0.0 s1];
   [eta,f] = ode45(@fcnjh,etaspan,f0,options);
   phi1 = f(end,1);
   ds = -phi1*(s1 - s0)/(phi1 - phi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   phi0 = phi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
% -------- Plotting f and f' as functions of eta
FS = 'FontSize';
clf
plot(eta,f(:,1),'k',eta,f(:,2),'k-.')
grid on
st1 = 'Jeffrey-Hamel flow : ';
st2 = sprintf('Re = %3.1f , \\alpha = %3.1f',Re,alphag);
st = [st1 st2];
xlabel('\eta',FS,14,'FontWeight','Bold')
ylabel('f'' , f"',FS,14)
title(st,FS,14)
legend('f','f''')
% ------- Compute a table of f, f' and f'' 
etaspan = [0:0.05:1.0];
f0 = [1.0 0.0 s];
[eta,f] = ode45(@fcnjh,etaspan,f0,options);
fprintf('\n         eta        f          f''        f"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[eta f]');