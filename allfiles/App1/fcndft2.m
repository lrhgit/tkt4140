function y = fcndft2(x,a0,a,b)
% Beregn et trigonometrisk polynom for
% vektoren x når Fourier-koeffisientene er gitt.
% Disse koeffisientene kan beregnes av funksjonen fcndft
%
%        --- Inn ---
% x : Gitt vektor med ekvidistante abscisser
%      på intervallet [0, 2*pi] (vilkårlig antall)
% a0 : 1. ledd i cosinusrekke fra fcndft
% a  : vektor med resten av koeffisientene for cosinusleddene
%      fra fcndft
% b:  vektor med koeffisientene for sinusleddene fra fcndft
%        --- Ut ---
% y: vektor med de interpolerte verdiene for vektoren x
%
%            === Referanse =
%     Uri M. Asher & Chen Greif :
%   " A First Course in NUMERICAL METHODS ",
%     SIAM 2011, kapittel 13 
%
l = length(a);
y = (a0/2)*ones(size(x));
for k = 1:l-1
    y = y + a(k)*cos(k*x) + b(k)*sin(k*x);
end
y = y + (a(l)/2)*cos(l*x);   
