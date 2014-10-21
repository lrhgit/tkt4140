% program timo2d
% Timoshenkos versjon.Dimensjonsløse størrelser
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen varierer lineært.
% Tilnærmet løsning der vi neglisjerer
% konstantene C(3) og C(4);
clear
% Lengder i inches og krefter i lbf.
R = 360; % Tankradius 
d = 312; % Høyde
x0 = 104; % Overhøyde
t1 = 3.5; % veggtykkelse øverst
t2 = 14; % veggtykkelse nederst
ny = 0.25;  % Poissons tall
ga = 0.03613; % Egenvekt [lbf/inch^3]
E = 30.0e6; % E-modul [lbf/inch^2]
a = t1/x0;
ro = (12*(1 - ny^2)*(x0/(a*R))^2)^0.25;
x1 = (x0 + d)/x0;
z1 = 2*ro*sqrt(x1);
K = ga*R^2/(E*a);
sx1 = sqrt(x1);

[berz1,beiz1] = berbei(z1);
[berdz1,beidz1] = berdbeid(z1);

J = z1*(berz1*berdz1 + beiz1*beidz1);

C1 = -2*(sx1*beidz1 - ro*d*berz1/x0)/J;
C2 = -2*(sx1*berdz1 + ro*d*beiz1/x0)/J;
% [kerz2,keiz2] = kerkei(z2);
% [kerdz2,keidz2] = kerdkeid(z2);

% A(1,1) = berdz2; A(1,2) = beidz2;
% A(1,3) = kerdz2; A(1,4) = keidz2;

% A(2,1) = -(2*berdz2 + z2*beiz2);
% A(2,2) = (-2*beidz2 + z2*berz2); 
% A(2,3) = -(2*kerdz2 + z2*keiz2);
% A(2,4) = (-2*keidz2 + z2*kerz2);

% y1 =H0/R; z1 = 2*ro*sqrt(y1);
% [berz1,beiz1] = berbei(z1);
% [berdz1,beidz1] = berdbeid(z1);
% [kerz1,keiz1] = kerkei(z1);
% [kerdz1,keidz1] = kerdkeid(z1);

% A(3,1) = 4*berdz1 - 0.5*z1^2*beidz1 + 2*z1*beiz1;
% A(3,2) = 4*beidz1 + 0.5*z1^2*berdz1 - 2*z1*berz1;
% A(3,3) = 4*kerdz1 - 0.5*z1^2*keidz1 + 2*z1*keiz1;
% A(3,4) = 4*keidz1 + 0.5*z1^2*kerdz1 - 2*z1*kerz1;

% A(4,1) = 2*beidz1 - z1*berz1 ;
% A(4,2) = -(2*berdz1 + z1*beiz1);
% A(4,3) = 2*keidz1 - z1*kerz1;
% A(4,4) = -(2*kerdz1 + z1*keiz1);

dy = 0.1; kmax = 1 + round(1/dy);
w = zeros(kmax,1); y = w; Mx = w; dw = w;
% dw = w;
for k = 1:kmax
    y(k) = (k - 1)*dy;
    x = d*(1- y(k))/x0 + 1;
    z = 2*ro*sqrt(x);
    [berz,beiz] = berbei(z);
    [berdz,beidz] = berdbeid(z); 
%     [kerz,keiz] = kerkei(z);
%     [kerdz,keidz] = kerdkeid(z);
    wp = - (x - 1)/x;
    w(k) = K*((C1*berdz - C2*beidz)/sqrt(x)+ wp);
    
    u1 = C1*(2*berdz + z*beiz);
    u2 = C2*(-2*beidz + z*berz);
    % u3 = -C(3)*(2*kerdz + z*keiz); u4 = C(4)*(-2*keidz + z*kerz);
    dw(k) = K*(-0.5*(u1 + u2 )/x^1.5 - 1/x^2)/x0;
    
    K2 = ga*a^2*R^2*x0/(12*(1 - ny^2));
    v1 = C1*(-z^2*beidz + 4*z*beiz+ 8*berdz);
    v2 = -C2*(z^2*berdz - 4*z*berz + 8*beidz);
    Mx(k) = -K2*(sqrt(x)*(v1 + v2)*0.25 + 2);
end
