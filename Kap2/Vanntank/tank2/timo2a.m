% program timo2a
% Timoshenkos versjon med asymptotiske uttrykk
% for Kelvinfunksjonene.
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen varierer lineært.
% Tilnærmet løsning der vi neglisjerer
% konstantene C(3) og C(4);
clear
% Lengder i inches og krefter i lbf.
R = 360; % Tankradius 
d = 312; % Høyde
x0 = 104; % Overhøyde
t1 = 3.5; % veggtykkelse øvers
t2 = 14; % veggtykkelse nederst
ny = 0.25;  % Poissons tall
ga = 0.03613; % Egenvekt [lbf/inch^3]
E = 30.0e6; % E-modul [lbf/inch^2]
a = t1/x0;
ro = (12*(1 - ny^2)/(a*R)^2)^0.25;
x1 = x0 + d;
z1 = 2*ro*sqrt(x1);
K = ga*R^2/(E*a);
N = sqrt(2*pi*z1)*exp(-z1/sqrt(2));
arg1 = z1/sqrt(2) + pi/8; arg2 = z1/sqrt(2) - pi/8; 
sx1 = sqrt(x1);

% C1a og C2a er konstanter ved bruk av asymptotiske
% uttrykk for Kelvinfunksjonene.
C1a = K*sqrt(2)*N*(-sx1*sin(arg1) + ro*d*cos(arg2))/(ro*sx1);
C2a = - K*sqrt(2)*N*(sx1*cos(arg1) + ro*d*sin(arg2))/(ro*sx1);

% === Beregning basert på asymptotiske uttrykk ===

dy = 0.05; kmax = 1 + round(1/dy);
w = zeros(kmax,1); y = w; Mx = w; dw = w;
for k = 1:kmax
    y(k) = (k - 1)*dy;
    x = d*(1- y(k)) + x0;
    z = 2*ro*sqrt(x);     
    wp = - K*(x - x0)/x;
    arg1 = z/sqrt(2) + pi/8; arg2 = z/sqrt(2) - pi/8;
    fac = exp(z/sqrt(2))/sqrt(2*pi*z);
    w(k) = fac*(C1a*cos(arg1) - C2a*sin(arg1))/sqrt(x)+ wp;
    
    u1 = C1a*(2*cos(arg1) + z*sin(arg2));
    u2 = C2a*(-2*sin(arg1) + z*cos(arg2));
    dw(k) = -0.5*fac*(u1 + u2 )/x^1.5 - K*x0/x^2;

    K2 = E*a^3/(48*(1 - ny^2));
    v1 = C1a*(-z^2*sin(arg1) + 4*z*sin(arg2)+ 8*cos(arg1));
    v2 = -C2a*(z^2*cos(arg1) - 4*z*cos(arg2)+ 8*sin(arg1));
    Mx(k) = -K2*sqrt(x)*fac*(v1 +v2) - ga*a^2*R^2*x0/(6*(1-ny^2));
    
end
