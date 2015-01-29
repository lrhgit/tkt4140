% Program plot21
% Finner startverdier for ligning 2.1
% i avsnitt 2.2 ved å plotte phi(s)
% som funksjon av s.
% Ser mer generelt på ligningen
%     y''(x)= a*y(x)^2 , a > 0
% En av løsningene er y = 6/(a*(1+x)^2)
% med y(0) = 6/a , y(1) = 3/(2*a)
% Denne løsningen har y'(0) = -12/a
% Beregner y(0) som funksjon av s
% der s = y(0)
%
clear; clear global a;
global a;
a = 2.5;
y0 = 6/a;
y1 = 3/(2*a);
sstart = -30; send = -4.5; antall = 10;
s = linspace(sstart,send,antall);
phi = s;
options= odeset('RelTol', 1.0e-5);
sspan = [0 1.0];
for n = 1:antall
   ys = [y0 s(n)];
   [x,y] = ode45(@fcn21,sspan,ys,options);
   phi(n)=y(end,1)-y1 ;
end
clf
FS = 'FontSize'; FW = 'FontWeight';
plot(s,phi)
grid
xlabel('s',FS,14)
ylabel('\phi','Rotation',0,FS,14,FW,'Bold')
st = sprintf('Plott av \\phi : a = %5.2f \n',a);
title(st,FS,14,FW,'Bold')
shg