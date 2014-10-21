% === tretrap===
% Beregner den analytiske løsningen
% av varmeledningsproblemet i det sammensatt
% trekantprofilet i avsnitt 3.2.2 i kompendiet.
clear
% === Data ===
La = 0.1; Da = 0.01; ha = 80; ka = 160; % trekantdelen
Lb = 0.1; Db = 0.02; db = Da; hb = 80; kb = 40; % trapesdelen

% === Koeffisienter ===
bta = sqrt(2*ha*La^2/(Da*ka));
btb = sqrt(2*hb*Lb^2/(Db*kb));
a = db/Db;
f = (bta/btb)*(ka/kb)*sqrt(a);

z0 = sqrt(a)*2*btb/(1 - a); z1 = 2*btb/(1-a); z1a = 2*bta;
K1z0 = besselk(1,z0); I0z1 = besseli(0,z1); 
K0z1 = besselk(0,z1); I1z0 = besseli(1,z0); 
K0z0 = besselk(0,z0); I0z0 = besseli(0,z0);
I0z1a = besseli(0,z1a); I1z1a = besseli(1,z1a);
J = I0z1a*(I0z1*K1z0+ I1z0*K0z1) + f*I1z1a*(I0z1*K0z0 - I0z0*K0z1);
A = (K0z0*I1z0 + I0z0*K1z0)/J;
B = (I0z1a*K1z0 + f*I1z1a*K0z0)/J;
C = (I0z1a*I1z0 - f*I1z1a*I0z0)/J;

% === Trekant ===
% For x = 0 : dtheta/dx = bta*bta*theta
% Dersom temperaturgradienten her skal være lik den 
% som vi får i den numeriske løsningen,
% må vi derfor bruke verdien av bta fra den numeriske løsningen
% i dtheta.
btan2 = 2.0;
h = 0.1;
x = [0 : h : 1.0]';
z = sqrt(x)*2*bta;
theta1 = A*besseli(0,z);
kmax = length(x);
dtheta1 = zeros(kmax,1);
dtheta1(1) = btan2*A;
for k = 2: kmax
    dtheta1(k) = A*btan2*besseli(1,z(k))./sqrt(x(k));
end
y = x*0.5;
fprintf('   x          theta             theta'' \n\n');
fprintf('%6.2f %18.12f %18.12f \n',[y,theta1,dtheta1]');

% === Trapes ===
% Må også her multiplisere med 2 for å få rett løsning
% for temperaturgradienten.
h = 0.1;
x = [0 : h : 1.0]';
t1 = sqrt(a +(1 - a)*x);
z = t1*2*btb/(1-a);
theta2 = B*besseli(0,z) + C*besselk(0,z);
dtheta2 = 2*btb*btb*(B*besseli(1,z) - C*besselk(1,z))./((1-a)*z);
y = x*0.5 + 0.5;
fprintf('%6.2f %18.12f %18.12f \n',[y,theta2 2*dtheta2]');


 