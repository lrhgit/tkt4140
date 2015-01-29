% ==== program lwdifferr ====
%
% Computes the diffusion error of the Lax-Wendroff 
% scheme for the advection equation
% as a function of the Courant number and the 
% phase-angle delta
% 
clear; close;
dvecg = linspace(0,180);
drad = dvecg*pi/360;% delta/2
sdel4 = sin(drad).^4;
ylim([0 1]);
LW = 'LineWidth';
hold on
for C = [0.25 0.5 0.8]
    C2 = C*C;
    G = sqrt(1 -4*C2*(1 - C2)*sdel4);
    plot(dvecg,G,LW,1);
end
hold off
grid
title('Diffusjonsfeil \epsilon{_D} for Lax-Wendroff skjema');
xlabel('\delta (grader)');


    