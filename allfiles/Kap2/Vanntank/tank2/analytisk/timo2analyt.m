% program timo2analyt
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen varierer lineært
% Lager tabeller for fysiske verdier
% basert på data nedenfor.
% Data er tatt fra
% Timoshenkos "Plates % Shells",
% avsnitt 118.
% 
clear
R = 360*0.0254; % Tankradius(m)
H = 312*0.0254;  % Høyde(m)
t0 = 14*0.0254;  % Veggtykkelse bunn (m)
t1 = 3.5*0.0254;  % Veggtykkelse topp
ny = 0.25;  % Poissons tall
E = 2.0e10;  % E-modul (Pa) . Ikke gitt hos Timoshenko
ga = 9807.4;  % Egenvekt (N/m^3)

a = (t0 - t1)/t0;
b4 = 3*(1 - ny^2)*H^4/(R*t0)^2;
b = (3*(1 - ny^2)*H^4/(R*t0)^2)^0.25;
ro = (b/a)*sqrt(2);
z0 = 2*ro;
K1 = ga*H*R^2/(E*t0); K2 = 0.25*ga*H^3/b4;
A = zeros(4,4); B = zeros(4,1); C = B;

[berz0,beiz0] = berbei(z0);
[berdz0,beidz0] = berdbeid(z0);
[kerz0,keiz0] = kerkei(z0);
[kerdz0,keidz0] = kerdkeid(z0);

A(1,1) = berdz0; A(1,2) = beidz0;
A(1,3) = kerdz0; A(1,4) = keidz0;

A(2,1) = -(2*berdz0 + z0*beiz0);
A(2,2) = (-2*beidz0 + z0*berz0); 
A(2,3) = -(2*kerdz0 + z0*keiz0);
A(2,4) = (-2*keidz0 + z0*kerz0);

y1 = t1/t0; z1 = z0*sqrt(y1);
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

fac = (a - 1)/a;
B(1) = 1;
B(2) = -2*fac;
B(3) = 4*fac/sqrt(y1);
B(4) = 0;
C = A\B;
dx = 0.1; kmax = round(1/dx);
w = zeros(kmax,1);
dw = w; m = w; v = w; x = w; 
for k = 1:kmax + 1
    x(k) = (k - 1)*dx;
    y = 1 - a*x(k);
    z = 2*ro*sqrt(y); 
    [berz,beiz] = berbei(z);
    [berdz,beidz] = berdbeid(z); 
    [kerz,keiz] = kerkei(z);
    [kerdz,keidz] = kerdkeid(z);
    wp = -(a + y - 1)/(a*y);
    w(k) =K1* ((C(1)*berdz + C(2)*beidz + C(3)*kerdz + C(4)*keidz)/sqrt(y)+ wp);
    u1 = -C(1)*(2*berdz + z*beiz); u2 = C(2)*(-2*beidz + z*berz);
    u3 = -C(3)*(2*kerdz + z*keiz); u4 = C(4)*(-2*keidz + z*kerz);
    dw(k) = -(K1*a/H)*(0.5*(u1 + u2 + u3 + u4)/y^1.5 + fac/y^2);
    u1 = C(1)*(4*berdz -0.5*z^2*beidz + 2*z*beiz);
    u2 = C(2)*(4*beidz + 0.5*z^2*berdz - 2*z*berz);
    u3 = C(3)*(4*kerdz -0.5*z^2*keidz + 2*z*keiz);
    u4 = C(4)*(4*keidz + 0.5*z^2*kerdz - 2*z*kerz);
    m(k) = -K2*(0.5*sqrt(y)*(u1 + u2 + u3 + u4) - 2*fac)*a^2;
    u1 = C(1)*(2*beidz - z*berz);
    u2 = -C(2)*(2*berdz + z*beiz);
    u3 = C(3)*(2*keidz - z*kerz);
    u4 = -C(4)*(2*kerdz + z*keiz);
    v(k) = (K2/H)*a^3*0.5*ro^2*sqrt(y)*(u1 + u2 + u3 + u4);
end
s1 = '    x         W            W''';
s2 = '               M            V \n \n';
fprintf([s1,s2]);
fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw m v]');