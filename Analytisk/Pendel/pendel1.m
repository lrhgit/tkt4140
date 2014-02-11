% Program pendel1
% Analytisk løsning av enkel pendel
clear
thet0g = input('Les theta0 i grader: ');
n = input('Les antall pkt. n: ');
L = input('Les pendel-lengde L i meter: ');

thet0 = thet0g*pi/180;% theta(0) i radianer
% T = fjerdeparten av perioden ( dimensjonsløs tid)
% for store utslag.

k = sin(thet0/2);
k2 = k*k;
T = ellipke(k2);
t = linspace(0,T,n);
tm = t(end:-1:1);

%Løsning for ikke-lineær ligning
[sn,cn]=ellipj(tm,k2);
theta = 2*asin(k*sn); 
dtheta = -2*k*cn; 
% --- Tabell. Fysiske verdier ---
g = 9.81; fac = sqrt(L/g);
theta = theta*180/pi;
t = fac*t; % Fysisk tid
dtheta = dtheta/fac; % rad/s
fprintf('\n tid(s)    theta(grd)    dtheta(rad/s)\n\n');
for l = 1:length(t)
    fprintf('%7.4f   %12.5e   %12.5e \n',t(l),theta(l),dtheta(l));
end;
% Plotting av resultater
plot(t,theta)
grid on
FS = 'FontSize';
st = sprintf('\\theta_0 = %4.2f\\circ',thet0g);
title(st,FS,13);
xlabel('t',FS,14);
ylabel('\theta','Rotation',0,FS,13);
shg