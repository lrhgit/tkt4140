%================ program analytorg ===================
% 
% Analytisk løsning av temperaturfordeling T(r,t) i et
% en kule med radius r der 0 <= r <= b.
% Kula med starttemperatur Tk blir sluppet i vann
% med konstant temperatur Tv.
% Denne løsningen forutsetter av H*b/K = 1
% der H er varmeovergangstallet og 
% K er varmeledningstallet.
%
clear
alfa = 0.04; % Termisk diffusivitet (cm^2/s)
t = 10; % Tid i sekund
b = 5; % radius av kule
r = 2.5; % radius i cm
Tk = 300.0;% Starttemperatur i kula
Tv = 20; % Temperatur i vannet.(konstant)
mmax = 201; % Antall ledd
fprintf('Antall ledd : %8.0f\n',mmax);
fac1 = pi*0.5/b;
tegn = 1;
s = 0;
for m = 1:2:mmax 
   bm = m*fac1;
   bmr = r*bm;
   temp = exp(-alfa*bm^2*t)/m^2;
   ledd = tegn*sin(bmr)*temp;
   s = s + ledd;
   tegn = -tegn;
end
T = Tv + 8*b*(Tk - Tv)*s/(pi^2*r)
% r = (0:dr:5)';
% fprintf(' Fouriertall : %8.3f\n',D);
% fprintf(' Tidsinkrement(sekund) : %8.3e\n',dtau);
% fprintf(' Antall tidskritt : %6d\n',ntot);
% fprintf(' Forløpet tid(sekund) : %8.3f\n',tau);
% fprintf('\n    r         T\n\n');
% fprintf('   %4.2f  %10.4f\n',[r Tnew]');
% plot(r,Tnew)
% grid
% shg