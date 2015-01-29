%========== Program kule3 ==========
% 
% Beregner temperaturfordeling i et
% en kule.
% Bruker FTCS-skjemaet.
% 2.ordens bakoverdifferanser
% i randbetingelsen for r = R = 5cm
% Bruker falskt punkt  for r = 0.
% I dette tilfellet er D gitt sammen med antall skritt.
%
clear
b = 5;
m = 4; % Antall deler
alfa = 0.04; % Termisk diffusivitet
dr = b/m;   % r-skritt i cm
D = 0.3; % Numerisk Fourier-tall
dtau = D*dr^2/alfa;   % Tidsinkrement i sekund
tauend = 700; % Stopptid
ntot = round(tauend/dtau); % antall tidskritt
K = 0.1; % Varmeledningstall
H = 0.02 ; %Varmeovergangstall
dra = dr*H/K;
Tk = 300;
Told = Tk*ones(m+1,1);% Initialverdier
Tnew = Told;
Tv = 20; % temperatur i vannet.
r = (0:dr:5)';
%Told(m+1) = (4*Told(m) -Told(m-1)+ 2*dra*Tv)/(3 + 2*dra); % 2.ordens 
for k = 1:ntot 
   for j = 2 :m
      Tnew(j)= (1 - 2*D)*Told(j)+ D*((1-1/(j-1))*Told(j-1)+ (1+1/(j-1))*Told(j+1)); 
   end
   Tnew(1) = (1 - 6*D)*Told(1) + 6*D*Told(2);
   if abs(Tnew(1)) > 1.1*Tk
      break;
   end
   Tnew(m+1) = (4*Tnew(m) -Tnew(m-1)+ 2*dra*Tv)/(3 + 2*dra); % 2.ordens
   Told = Tnew;
   tau = k*dtau;
%     fprintf(' Forløpet tid(timer) : %6d\n',tau);
%     fprintf('   %4.2f  %10.4f\n',[x Tnew]');
end
fprintf(' Fouriertall : %8.3f\n',D);
fprintf(' r-inkrement : %10.3e\n',dr);
fprintf(' Tidsinkrement(sekund) : %8.3e\n',dtau);
fprintf(' Antall tidskritt : %6d\n',ntot);
fprintf(' Forløpet tid(sekund) : %8.3f\n',tau);
% fprintf('\n    r         T\n\n');
% fprintf('   %4.2f  %10.4f\n',[r Tnew]');
plot(r,Tnew)
grid
shg