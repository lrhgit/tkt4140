% ================== Program vtrapes ==================
% Computes the analytical solution
% of the heat-conduction problem for the
% trapezoidal profile in section 3.2 in the compendium
% 
clear
z0 = 4.0*sqrt(2.0); z1 = 8.0; g = sqrt(2.0)/8.0;
K1z0 = besselk(1,z0); I0z1 = besseli(0,z1); 
K0z1 = besselk(0,z1); I1z0 = besseli(1,z0);  
K0z0 = besselk(0,z0); I0z0 = besseli(0,z0);
J = K1z0*I0z1 + K0z1*I1z0 + g*(K0z0*I0z1 - K0z1*I0z0);
A = (K1z0 + g*K0z0)/J;
B = (I1z0 - g*I0z0)/J;
x = (0 : 0.1 : 1.0)';
z = sqrt(32.0*(1+x));
theta = A*besseli(0,z) + B*besselk(0,z);
dtheta = 16.0*(A*besseli(1,z) - B*besselk(1,z))./z;
fprintf('   x       theta       theta'' \n\n');
fprintf('%6.2f %10.5f  %10.5f \n',[x,theta,dtheta]');
