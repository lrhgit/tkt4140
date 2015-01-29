% === Program plot22 ===
% Finner startverdier for ligning 2.1
% i avsnitt 2.2 ved å plotte phi(s)
% som funksjon av s.
% Ligningen er gitt ved :
%     y''(x)= (3/2)*y(x)^2 
% Beregner y(0) som funksjon av s
% der s = y(0)
%
clear;close;
y0 = 4.0;
sstart = -12; send = -6; antall = 20;
s = linspace(sstart,send,antall);
phi = s; % Allokerer plass
options= odeset('RelTol', 1.0e-5);
sspan = [0 1.0];
for n = 1:antall
   ys = [y0 s(n)];
   [x,y] = ode45(@fcn22,sspan,ys,options);
   phi(n)=y(end,1)- 1 ;
end
FS = 'FontSize'; FW = 'FontWeight';
plot(s,phi)
grid
xlabel('s',FS,14)
ylabel('\phi','Rotation',0,FS,14,FW,'Bold')
st = sprintf('Plott av \\phi \n');
title(st,FS,14,FW,'Bold')