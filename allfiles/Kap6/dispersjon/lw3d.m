% program lw3d
% 
clear; close;
cv = linspace(0,1,40);
dv = linspace(0,pi,40);
[x,y] = meshgrid(dv,cv);
svart = zeros(64,3);
colormap(svart);
z = sqrt(1-4*y.^2.*(1-y.^2).*sin(x*0.5).^4);
mesh(x,y,z);
xlim([0,pi]);
box on
grid on
title('Amplitude |G| som funksjon av C og \delta');
ylabel('C');
xlabel('\delta(rad.)');
