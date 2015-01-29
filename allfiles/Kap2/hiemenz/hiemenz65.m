% ============= Program hiemenz ======================
% Løsning av Falkner-Skan ligningen for beta = 1 (Hiemenz)
% ved bruk av skyteteknikk. (Se eksempel 2.5 i kompendiet)
% Hiemenz ligning er gitt ved:
%    f'''(eta) + f(eta)*f''(eta)+ (1- (f')^2) = 0
%    f(0) = 0, f'(0) = 0, f'(etainf)=1
% I programmet setter vi f(1)= f, f(2) = f'
% og f(3) = f'' .
% 
% Programmet leser en verdi for s0, s1 og etainf
%  s = f''(0) = f(3)(0)
%
clear
etainf = input(' etainf = ?');
s0 = input(' s0 = ?');
s1 = input(' s1 = ?');

fprintf('        etainf = %8.2f\n',etainf);
fprintf('        s0 = %8.5f s1 = %8.5f   \n\n',s0,s1);
eta0 = 0; 
espan = [eta0 etainf];
% Beregner fi0
f0 = [0.0 0.0  s0];
[eta,f] = ode45('fcnhiemenz',espan,f0);
fi0 = f(end,2) -1;
% Startverdier for iterasjonen
itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
options = odeset('RelTol',1.0e-5);
% Skriv overskrift for tabell
fprintf('        itr.      s          ds\n\n');
% Starter iterasjon
while(abs(ds) > epsi) && (it < itmax)
   it = it + 1;
   f0 = [0.0 0.0 s1];
      [eta,f] = ode45('fcnhiemenz',espan,f0,options);
   fi1 = f(end,2) -1;
   ds = -fi1*(s1 - s0)/(fi1 - fi0);
   s = s1 + ds;
   s0 = s1;
   s1 = s;
   fi0 = fi1;
   fprintf('%10d %12.6f %12.3e\n',it,s,ds);
end
% Beregner en tabell for f, f' og f'' 
deta = 0.25;
espan = (eta0:deta:etainf);
if(mod(etainf,deta)~=0)
   espan = [espan etainf];
end
f0 = [0.0 0.0 s];
[eta,f] = ode45('fcnhiemenz',espan,f0,options);
fprintf('\n         eta        f          f''        f"\n\n');
fprintf(' %12.2f %10.6f %10.6f % 13.5e\n',[eta f]');
% Plotter f' og f" som funksjon av eta
espan = [eta0 etainf];
f0 = [0.0 0.0 s];
[eta,f] = ode45('fcnhiemenz',espan,f0,options);
clf
plot(f(:,2),eta)
hold on
plot(f(:,3),eta,'-.')
grid on
ylabel('\eta','FontSize',14,'FontWeight','Bold','Rotation',0)
xlabel('f'' , f"','Fontsize',14)
title('Løsning av Hiemenz ligning','Fontsize',14)
legend('f''','f"')
hold off
shg