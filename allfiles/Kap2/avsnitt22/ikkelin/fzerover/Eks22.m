% =================== eks22 ===========================
% Programmet løser randverdi-problemet i 
% avsnitt 2.2 i kompendiet ved bruk av skyteteknikk.
% Bruker FZERO for nullpunktbestemmelse
% Ligning : y''(x)= (3/2)*y(x)^2
%           y(0) = 4, y(1) = 1
% =======================================================
clear
xspan = [0.0 1.0];
% Angir et intervall s0 som vi vet inneholder
% et nullpunkt for phi. Intervallet kan finnes ved å plotte
% phi(s). Dersom vi ikke kjenner et slikt intervall, 
% kan vi angi en startverdi i nærheten av et nullpunkt.
s0 = [-6.0 -9.0];
% --- Setter relativ nøyaktighet for s = 1.0e-5
options = optimset('Display','off','TolX',1.0e-5);
tic
s = fzero(@fcnphi,s0,options);
toc
fprintf('Beregnet verdi av y''(0) = %12.5e \n',s);
options = odeset('RelTol',1.0e-5);
% ---- Beregner en tabell for y og y'
xspan = [0:0.1:1.0];
y0 = [4.0 s];
options = odeset('RelTol',1.0e-5);
tic
[x,y] = ode45(@fcn22,xspan,y0,options);
toc
fprintf('\n           x         y           y''\n\n');
fprintf(' %12.2f %12.6f %12.6f\n',[x y]');
