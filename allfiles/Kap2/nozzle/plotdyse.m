% ============== program plotdyse =========================
% Finner startverdier for f'(0) i plan dyse.
% Beregner phi som funksjon av s der s = f'(0)
% Kaller fcndyse
% =======================================================
clear
sstart = 0; send = 2.5; numbers = 50;
s = linspace(sstart,send,numbers);
etainf = 4.0;
K = 1;
phi = s;
options= odeset('RelTol', 1.0e-5);
etaspan = [0 etainf];
for n = 1:numbers
   f0 = [0 ; s(n); 0];
   [eta,f] = ode45(@fcndyse,etaspan,f0,options);
   phi(n) = f(end,1) - K ;
end
%st = sprintf('Exercise 5.1 :   p{_1}/p{_0} = % 5.2f', p1p0);
clf
plot(s,phi)
grid
FS = 'FontSize'; FW = 'FontWeight';
xlabel('s',FS,14)
ylabel('\phi','Rotation',0,FS,14,FW,'Bold')
%title(st,FS,13)
shg