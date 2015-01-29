%========== PROGRAM kule ========
%    --- Eksempel 5.2 ---
% Beregner temperaturfordeling i et
% en kule med radius b som avkjøles i vann
% Starttemperatur = Tk
% Temperatur i vann  = Tv = konstant
% Beregning ved bruk av FTCS-skjemaet.
% r-skritt dr og det numeriske Fourier-tallet D
% gis og tidskrittet dtau beregnes.
% Stopptiden tauend gis.
%
% Bruker 2.ordens bakoverdifferanser
% i randbetingelsen for r = b = 5cm
% rb = 1 : Bruker randbetingelse med separat ligning
%          for r = 0. Stabil for D <1/3
% rb = 2 : Bruker 2.ordens foroverdifferanser  for r = 0.
%          Stabil for D < 1/2
%     
% Den analytiske løsningen beregnes av
% funksjonen kanalyt. Betingelse for bruk av
% den analytiske løsningen : H*b/K = 1
% der H = varmeovergangstallet og
% K = varmeledningstallet
%
clear
clear global alfa b Tk Tv;
global alfa b Tk Tv;
b = 5; % radius av kule
m = 50; % Antall deler
alfa = 0.04; % Termisk diffusivitet (cm^2/s)
dr = b/m;     % r-skritt i cm
D = 0.40; % Numerisk Fourier-tall
dtau = D*dr^2/alfa; % Tidsinkrement i sekund
tauend = 600;     % Stopptid i sekunder
ntot = round(tauend/dtau); % antall tidskritt
K = 0.1; % Varmeledningstall (W/cm/C)
H = 0.02 ; %Varmeovergangstall (W/cm^2/C)
dra = dr*H/K;
Ta = zeros(m+1,1); % Analytiske verdier
Tk = 300; % Starttemperatur i kula
Tv = 20; % Temperatur i vannet.
Told = Tk*ones(m+1,1);% Initialverdier
Tnew = Told;
rb = 2;
for k = 1:ntot 
   for j = 2 :m
      Tnew(j)= (1 - 2*D)*Told(j)+ D*((1-1/(j-1))*Told(j-1) + ...
               (1+1/(j-1))*Told(j+1)); 
   end
   if rb == 2
     Tnew(1) = (4*Tnew(2) - Tnew(3))/3; % 2. ordens
   elseif rb == 1
     Tnew(1) = (1 - 6*D)*Told(1) + 6*D*Told(2);
   end  
   Tnew(m+1) = (4*Tnew(m) -Tnew(m-1)+ 2*dra*Tv)/(3 + 2*dra); % 2.ordens
   Told = Tnew;
   tau = k*dtau;
%     fprintf(' Forløpet tid(timer) : %6d\n',tau);
%     fprintf('   %4.2f  %10.4f\n',[x Tnew]');
end
r = (0:dr:5)';
% Analytisk løsning
for k = 1:length(r)
    x = r(k);
    Ta(k) = analyt(x,tau);
end
    
fprintf(' Fouriertall : %8.3f\n',D);
fprintf(' Tidsinkrement(sekund) : %8.3e\n',dtau);
fprintf(' Antall tidskritt : %6d\n',ntot);
fprintf(' Forløpet tid(sekund) : %8.3f\n',tau);
fprintf('\n    r         T            Ta \n\n');
fprintf('   %4.2f  %10.2f  %10.2f \n',[r Tnew Ta]');
% plot(r,Tnew)
% grid
% shg