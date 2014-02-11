% ==== program lwdispserr ====
%
% Computes the dispersion error of the Lax-Wendroff 
% scheme for the advection equation
% as a function of the Courant number and the 
% phase-angle delta
% 
clear; close;
dvecg = linspace(0,180);
drad = dvecg*pi/180;% delta
drad2 = drad/2;% delta/2
sdel = sin(drad); % sin(delta)
sdel2 = 2*sin(drad2).^2;% 2*sin(delta/2)^2
%ylim([0 1]);
LW = 'LineWidth';
hold on
for C = [0.25 0.5 0.8]
    temp = atan2(C*sdel,1 - C^2*sdel2);
    derr = temp./(C*drad);
    plot(dvecg,derr,LW,1);
end
hold off
grid
title('Dispersjonsfeil \epsilon{_\phi} for Lax-Wendroff skjema');
xlabel('\delta (grader)');

    