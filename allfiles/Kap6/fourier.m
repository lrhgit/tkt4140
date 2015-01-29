% program fourier
% Fourierrekke for funksjonen f(x) = 1 - x^2
% på intervallet -1 <= x <= 1
% Vektoren v gir antall ledd i Fourierutviklingen for hvert
% av tilfellene, der length(v) angir antall tilfeller.
% Hvert tilfelle lagres i matrisa y der hver linje inneholder
% ordinatene for tilhørende x-verdi.
clear; close;
x = linspace(-1,1);
v = [1  3 10];
m = length(v);
y = zeros(m,length(x));
s = zeros(size(x));
s = 2/3;
fac = (2/pi)^2;
for N = 1:m
    n = v(N);
    s = 2*ones(size(x))/3;
    tegn = 1;
    for k = 1:n   
        s = s + tegn*fac*cos(k*pi.*x)/k^2;    
        tegn = -tegn;
    end
    y(N,:) = s;
end
% plot(y(1:m,:)');
% grid
surf(y)
shading interp
axis off ij