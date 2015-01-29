% Program test
%
% Som advecus1, men kompakt form av iterasjonsløkke
% med bruk av indeksvektor.
%
% Oppstrømskjema for adveksjonsligningen
% Genererer en av figurene i fig. 6.5 i manus
% Leser inn Courant-tallet C og antall tidskritt nitr
% Bruker dx = 0.01
% Eksempel på brukbare verdipar :
% C = 1 , nitr = 10 
% C = 0.5 , nitr = 20
% C = 0.25 , nitr = 40
%
clear; close;
n = 11; % Antall punkt
dx = 1/(n-1);
u = zeros(n,1);
%u0 = u;
x = u; v = u;
x = 0 : dx : 1.0;
% Startverdier for u.
% u = 1.0 for x < 0.5 ,ellers u = 0
u = 1.0*(x < 0.5);
C = input('Courant-tall = ?');
nitrmax = 0.5/(C*dx);
nitr = input('Antall iterasjoner = ?');
if nitr > nitrmax
    fprintf ('Max. iterasjoner er %4.0f \n',nitrmax);
    nitr = input('Antall iterasjoner = ?');
end
xfront = nitr*C*dx + 0.5; % Front of exact solution
u0 = 1.0*(x < xfront);% uexact = 1 for x < xfront
J = [2 : n]; % Indeksvektor
for itr = 1: nitr
    v = u;  
    u(J) = (1 - C)*v(J) + C*v(J-1);
end
FS = 'FontSize';
%plot(x,u,'k',x,u0,'k:');
plot(x,u,'k');
ylim([-0.5 1.5]);
xlabel('x',FS,14);
ylabel('u(x,t)',FS,14);
st = sprintf('C = %4.2f',C);
title(st,FS,14);
%subplot(312),plot(x,u,'k',x,u0,'k:');
%subplot(313),plot(x,u,'k',x,u0,'k:');
%ylim([-2.0 2.0]);