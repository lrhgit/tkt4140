% ==== program ampfacv2 ====
%
% Computes the amplification factor G of the 
% FTCS-scheme for diffusion equation
% as a function of the diffusion number D and the 
% phase-angle delta
%
clear all; close all; clc;
FS = 20; LW=3; set(0,'DefaultLineLineWidth',LW,'DefaultAxesFontName','Helvetica','DefaultAxesFontSize',FS);
dvecg = linspace(0,180,50);
drad = dvecg*pi/180;% delta i radianer
sdel = sin(drad/2).^2;
ylim([-1 1]);
hold on
% Ga
for D = [0.25 0.5]
    Ga = exp(-D*drad.^2);
    plot(dvecg,Ga,'k-.');
end
% G
for D = [0.25 0.35 0.5]
    G = 1 -4*D*sdel;
    plot(dvecg,G,'k');
end
hold off
grid
title('Forsterkningsfaktorer for diffusjonligningen');
xlabel('\delta (grader)');
ylabel('G{_a} , G');
legend('G','G{_a}');


    