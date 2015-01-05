function y = fcndct2(x,a0,a)
% Beregner et trigonometrisk polynom for
% vektoren x når Fourier-koeffisientene for den diskrete
% cosinus-transformasjonen er er beregnet fra fcndct.
%
%        --- Inn ---
% x : Gitt vektor med ekvidistante abscisser
%      på intervallet [0, pi] (vilkårlig antall)
% a0 : 1. ledd i cosinusrekke fra fcndct
% a  : vektor med resten av koeffisientene for cosinusleddene
%      fra fcndct.
%        --- Ut ---
% y: vektor med de interpolerte verdiene for vektoren x
%
%            === Referanse =
%     Uri M. Asher & Chen Greif :
%   " A First Course in NUMERICAL METHODS ",
%     SIAM 2011, kapittel 13 
%
n = length(a) - 1;
y = (a0/2)*ones(size(x));
for k = 1:n
    y = y + a(k)*cos(k*x);
end 
