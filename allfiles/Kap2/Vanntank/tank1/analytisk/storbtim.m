% program storbtim
% Avsnitt 2.4.1
% Beregner innspenningsmomentet m0,M0 og
% skjærkrafta v0, V0 for en sirkulær 
% vanntank der veggtykkelsen er konstant.
% Bruker her den tilnærmede løsningen
% for store beta-verdier
% Bruker Timoshenkos data
clear
R = 360*0.0254; % Tankradius(m)
H = 312*0.0254;  % Høyde(m)
t = 14*0.0254;  % Veggtykkelse (m)
ny = 0.25;  % Poissons tall
ga = 9807.4;  % Egenvekt (N/m^3)

b4 = 3*(1 - ny^2)*H^4/(R*t)^2;
b = b4^0.25;
%K1 = ga*H*R^2/(E*t);
K2 = 0.25*ga*H^3/b4;

fprintf('beta = %12.5e \n',b);
m0 = 2*b*(b - 1); % Dimensjonsløs
v0 = -2*b^2*(2*b - 1); % Dimensjonsløs
M0 = K2*m0;  % dimensjon N
V0 = K2*v0/H; % dimensjon N/m
fprintf('m0 = %14.4f \n',m0);
fprintf('v0 = %14.4f \n',v0);
fprintf('M0 = %12.3f \n',M0);
fprintf('v0 = %12.3f \n',V0);

