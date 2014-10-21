% ==== program udifferr ====
%
% Computes the diffusion error of the 1. order 
% upwind scheme for the advection equation
% as a function of the Courant number and the 
% phase-angle delta
% 
clear; close;
dvecg = linspace(0,180);
drad = dvecg*pi/360;% delta/2
sdel = sin(drad).^2;
ylim([0 1]);
LW = 'LineWidth';
hold on
for C = [0.25 0.5 0.8]
    G = sqrt(1 -4*C*(1 - C)*sdel);
    plot(dvecg,G,LW,1);
end
hold off
grid
title('Diffusjonsfeil \epsilon{_D} for oppstrømskjema');
xlabel('\delta (grader)');


    