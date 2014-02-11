% ======================= Program hiemenz ======================
% Solution of the Falkner-Skan equation for beta = 1 (Hiemenz)
% using a shooting-technique. (See example 2.5 in the compendium)
% Hiemenz's equation :
%    f'''(eta) + f(eta)*f''(eta)+ (1- (f')^2) = 0
%    f(0) = 0, f'(0) = 0, f'(etainf) = 1
% In the program we put f(1)= f, f(2) = f' and f(3) = f'' .
%
% The program reads a value of s0, s1 and etainf
%
%  s = f''(0) = f(3)(0)
% Using ode45 with structure-syntax.
%
clear
etainf = input(' etainf = ?');
s0 = input(' s0 = ?');
s1 = input(' s1 = ?');

fprintf('        etainf = %8.2f\n',etainf);
fprintf('        s0 = %8.5f s1 = %8.5f   \n\n',s0,s1);
eta0 = 0; 
espan = [eta0 etainf];
% Compute fi0
f0 = [0.0 0.0  s0];
[eta,f] = ode45(@fcnhiemenz,espan,f0);
fi0 = f(end,2) -1;
% Initial values for the iteration
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
% Heading of table
fprintf('        itr.      s          ds\n\n');
% Start of iteration
while(abs(ds) > epsi) && (it < itmax)
   it = it + 1;
   f0 = [0.0 0.0 s1];
      [eta,f] = ode45(@fcnhiemenz,espan,f0,options);
   fi1 = f(end,2) -1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
% Compute a table of f, f' and f'' 
deta = 0.25;
espan = [eta0 etainf];
f0 = [0.0; 0.0; s];
sol = ode45(@fcnhiemenz,espan,f0,options);
eta = (eta0: deta: etainf);
if(mod(etainf,deta)~=0)
   eta = [eta etainf];
end
fvalue = deval(sol,eta); 
fprintf('\n         eta        f          f''        f"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[eta' fvalue']');
% Plotting f' og f" as a function of eta
eta = sol.x;
fvec = sol.y;
clf
plot(fvec(2,:)',eta,fvec(3,:)',eta,'-.')
grid on
ylabel('\eta','FontSize',14,'FontWeight','Bold','Rotation',0)
xlabel('f'' , f"','Fontsize',14)
title('Løsning av Hiemenz ligning','Fontsize',14)
legend('f''','f"')
shg