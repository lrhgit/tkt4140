function [a0,a] = fcndct(y)
% Diskret cosinus-transformasjon av m ekvidistante punkt
% (xj,yj), j = 0,1,.. m-1 der x(j) = pi(j + 1/2)/m
% og m er et partall.
%
%        --- Input ---
% y: vektor med  m y-verdier y0,y1,y2,..,yn
%         (kolonne eller linjevektor)
%        --- Output ---
% a0 : 1. ledd i cosinusrekke
% a  : vektor med koeffisienter for cosinusleddene.
%      Antall ledd er m/2 
% a er linjevektor.
%
%            === Referanse =
%     Uri M. Asher & Chen Greif :
%   " A First Course in NUMERICAL METHODS ",
%     SIAM 2011, kapittel 13 
%
y = y(:); % Konverter y til kolonnevektor.
m = length(y); n = m -1;
iv = 0: n;
x = pi*(iv + 1/2)/m;
a0 = 2*sum(y)/m;
% Løkke for beregning av a(k)
for k = 1: n 
    co = cos(k*x);
    a(k) = 2*(co*y)/m;
end