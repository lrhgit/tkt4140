% program timostorb
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen varierer lineært
% Beregner m(0) og v(0) for store beta-verdier
% (SI og engelsk) 
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
m0 = 2*b*(b-1) - a*(1-a);
v0 = -2*b^2*((2*b - 1) -(b-1)*a/b);
M0 = m0*K2;
V0 = v0*K2/H;
% s1 = '    x         W            W''';
% s2 = '               M            V \n \n';
% fprintf([s1,s2]);
% fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw m v]');