% ================= Program plotjh =====================
% Find starting-values for a root-finder in 
% for the Jeffrey-Hamel's equation, 
% by plotting phi(s) where s = f''(0).
% Using function fcnjh
% ======================================================
clear
global Re alpha;
alphag = 10;
alpha = alphag*pi/180;
Rea = 10.1797;
Re = Rea/alpha;
fprintf('Re*alfa = %15.7e \n',Rea);
fprintf('Re      = %15.7e \n',Re);
sstart= -7.5; send = -5.5 ; numbers = 30;
s = linspace(sstart,send,numbers);
phi = s;
options= odeset('RelTol', 1.0e-5);
etaspan = [0 1.0];
for n=1:numbers
   f0 = [1.0 0.0 s(n)];
   [eta,f] = ode45(@fcnjh,etaspan,f0,options);
   phi(n) = f(end,1);
end
FS = 'FontSize';
plot(s,phi)
grid on
st = sprintf('\\phi-plot: Re = %3.1f , \\alpha = %3.1f',Re,alphag);
xlabel('s',FS,14)
ylabel('\phi',FS,14,'FontWeight','Bold')
title(st,FS,14)
shg