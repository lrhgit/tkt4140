% ================= Program plotjh2 =====================
% Find starting-values for a root-finder in 
% for the Jeffrey-Hamel's equation for the 
% restricted case alpha^2 << Re*alpha,
% by plotting phi(s) where s = f''(0).
% Using function fcnjh2
% ======================================================
clear
global Rea;
Rea = -100.0;
fprintf('Re*alfa = %15.7e \n',Rea);
sstart= -0.0005; send = 0.0 ; numbers = 30;
s = linspace(sstart,send,numbers);
phi = s;
options= odeset('RelTol', 1.0e-5);
etaspan = [0 1.0];
for n=1:numbers
   f0 = [1.0 0.0 s(n)];
   [eta,f] = ode45(@fcnjh2,etaspan,f0,options);
   phi(n) = f(end,1);
end
FS = 'FontSize';
plot(s,phi)
grid on
st = sprintf('\\phi-plot: Re*alpha = %10.5f ',Rea);
xlabel('s',FS,14)
ylabel('\phi',FS,14,'FontWeight','Bold')
title(st,FS,14)
shg