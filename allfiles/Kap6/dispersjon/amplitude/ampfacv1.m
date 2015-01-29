% ==== program ampfac ====
%
% Computes the amplification factor G of the 
% FTCS-scheme for diffusion equation
% as a function of the diffusion number D and the 
% phase-angle delta
%
clear; close;
dvecg = linspace(0,180,50);
drad = dvecg*pi/180;% delta i radianer
sdel = sin(drad/2).^2;
ylim([-1 1]);
hold on
for D = [0.25 0.5]
    G = 1 -4*D*sdel;
    Ga = exp(-D*drad.^2);
    plot(dvecg,G,'k');
    plot(dvecg,Ga,'k-.');
end
hold off
grid
title('Forsterkningsfaktorer for diffusjonligningen');
xlabel('\delta (grader)');
ylabel('G{_a} , G');
legend('G','G{_a}');


    