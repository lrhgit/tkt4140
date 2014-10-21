% program ksum
% Verifiserer analytisk løsning av kule
% for t = 0.
% lam = r/b slik at 0 < lam < 1
% Summen skal bli : s = lam *pi^2/8
%
lam = 0.8;
fac = lam*pi/2;
s = 0;
tegn = 1;
mmax = 4001;
for m = 1:2:mmax
 arg = m*fac;
 ledd = sin(arg)/m^2;
 s = s + tegn*ledd;
 tegn = - tegn;
end
st = lam*pi^2/8;
fprintf( 'Teoretisk sum = %12.5e \n',st);
fprintf( 'Beregnet sum = %12.5e \n',s);
