% Program ex134
%
% Eksempel 13.4 i referansen nedenfor.
%
% --- Trigonometrisk interpolering ---
% Gitt m ekvidistante punkt (x1,y1),(x2,y2),..,(xn,yn)
% fra en funksjon y = f(x)
% Bruker funksjonen fcndft til å beregne koeffisientene
% a og b for det trigonometriske polynomet som
% går gjennom de gitte punktene.
% Interpoleringen utføres av funksjonen fcndft2
%
%            === Referanse =
%     Uri M. Asher & Chen Greif :
%   " A First Course in NUMERICAL METHODS ",
%     SIAM 2011, kapittel 13 
%
clear; close;
m  = 16; l = m/2;
x = 0: pi/l : (m-1)*pi/l; % abscisser på [0,2*pi]
j = 0: m-1;
t = 3*j/m - 1;
% y-verdier på intervallet [-1 , 2]
y = t.^2.*(t+1).^2.*(t-2).^2 - ...
    exp(-t.^2).*sin(t+1).^2.*sin(t-2).^2;
% Beregner koeffisientene a0,a og b;
[a0,a,b] = fcndft(y);
% Lager et fint nett på [0,2*pi] for interpolering.
xi = 0 : 0.01*pi : 2*pi;
% Beregner de interpolerte y-verdiene med bruk av dft2
yi = fcndft2(xi,a0,a,b);
tt = (3/(2*pi))*xi -1;
yex = tt.^2.*(tt+1).^2.*(tt-2).^2 - ...
    exp(-tt.^2).*sin(tt+1).^2.*sin(tt-2).^2;
plot(tt,yex,'b--',tt,yi)
axis([-1 2 -2 6])
