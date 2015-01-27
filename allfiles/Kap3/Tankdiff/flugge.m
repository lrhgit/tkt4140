% program flugge
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen varierer lineært
clear
R = 2.75; % Tankradius
H = 3.65; % Høyde
H0 = 1.35; % Overhøyde
t2 = 0.275; % veggtykkelse nederst
ny = 0.0;  % Poissons tall
a = t2/(H + H0);
ro = (12*(1 - ny^2)/a^2)^0.25;
y2 = (H + H0)/R;
z2 = 2*ro*sqrt(y2);
A = zeros(4,4); B = zeros(4,1); C = B;

[berz2,beiz2] = berbei(z2);
[berdz2,beidz2] = berdbeid(z2);
[kerz2,keiz2] = kerkei(z2);
[kerdz2,keidz2] = kerdkeid(z2);

A(1,1) = berdz2; A(1,2) = beidz2;
A(1,3) = kerdz2; A(1,4) = keidz2;

A(2,1) = -(2*berdz2 + z2*beiz2);
A(2,2) = (-2*beidz2 + z2*berz2); 
A(2,3) = -(2*kerdz2 + z2*keiz2);
A(2,4) = (-2*keidz2 + z2*kerz2);

y1 =H0/R; z1 = 2*ro*sqrt(y1);
[berz1,beiz1] = berbei(z1);
[berdz1,beidz1] = berdbeid(z1);
[kerz1,keiz1] = kerkei(z1);
[kerdz1,keidz1] = kerdkeid(z1);

A(3,1) = 4*berdz1 - 0.5*z1^2*beidz1 + 2*z1*beiz1;
A(3,2) = 4*beidz1 + 0.5*z1^2*berdz1 - 2*z1*berz1;
A(3,3) = 4*kerdz1 - 0.5*z1^2*keidz1 + 2*z1*keiz1;
A(3,4) = 4*keidz1 + 0.5*z1^2*kerdz1 - 2*z1*kerz1;

A(4,1) = 2*beidz1 - z1*berz1 ;
A(4,2) = -(2*berdz1 + z1*beiz1);
A(4,3) = 2*keidz1 - z1*kerz1;
A(4,4) = -(2*kerdz1 + z1*keiz1);

B(1) = sqrt(y2)*H/(H + H0);
B(2) = 2*H0/(R*sqrt(y2));
B(3) = - 4*H0/(R*sqrt(y1));
B(4) = 0;
C = A\B;
dx = 0.1; kmax = 1 + round(1/dx);
w = zeros(kmax,1); dw = w; x = w; x = w;
for k = 1:kmax
    x(k) = (k - 1)*dx;
    y = (H*(1- x(k)) + H0)/R;
    z = 2*ro*sqrt(y); 
    [berz,beiz] = berbei(z);
    [berdz,beidz] = berdbeid(z); 
    [kerz,keiz] = kerkei(z);
    [kerdz,keidz] = kerdkeid(z);
    wp = - (y - H0/R)/y;
    w(k) = (C(1)*berdz + C(2)*beidz + C(3)*kerdz + C(4)*keidz)/sqrt(y)+ wp;
    u1 = -C(1)*(2*berdz + z*beiz); u2 = C(2)*(-2*beidz + z*berz);
    u3 = -C(3)*(2*kerdz + z*keiz); u4 = C(4)*(-2*keidz + z*kerz);
    dw(k) = 0.5*(u1 + u2 + u3 + u4)/y^1.5 - H0/(R*y^2);
    
end