function [a0,a,b] = fcndft(y)
% Diskret Fourier transformasjon av m ekvidistante punkt
% (xj,yj), j = 0,1,.. m-1 der x(j) = (2*pi/m)*j
% og m er et partall.
%
%        --- Inn ---
% y: vektor med  m y-verdier y0,y1,y2,..,yn
%         (kolonne eller linjevektor)
%        --- Ut ---
% a0 : 1. ledd i cosinusrekke
% a  : vektor med koeffisienter for cosinusleddene.
%      Antall ledd er m/2
% b: vektor med koeffisienter for sinusleddene.
%     Antall ledd er m/2 - 1 
% a og b er linjevektorer
%
%            === Referanse =
%     Uri M. Asher & Chen Greif :
%   " A First Course in NUMERICAL METHODS ",
%     SIAM 2011, kapittel 13 
%
y = y(:); % Konverter y til kolonnevektor
m = length(y); l = m/2;
pi2m = 2*pi/m;
x = 0 : pi2m : (m -1)*pi2m;
a0 = sum(y)/l;
% Løkke for beregning av a(k) ok b(k)
for k = 1:l
    co = cos(k*x);
    si = sin(k*x);
    a(k) = (co*y)/l;
    if k < l
        b(k) = (si*y)/l;
    end
end