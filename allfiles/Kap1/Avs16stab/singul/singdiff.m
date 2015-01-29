% program singdiff
% plotter løsning av singulær diff.ligning
% i avsnitt 1.6.2
clear; clf;
x = linspace(0,1,300);
hold on
for epsilon = [0.03 0.1 0.2 0.4]
    fac = exp(1/epsilon) - 1;
    y = (exp(x/epsilon)-1)/fac;
    plot(x,y,'LineWidth',1.0);
end
grid on
xlabel('x','FontSize',14,'FontWeight','Bold')
ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
title('Grensesjikt ved x = 1','FontSize',14,'FontWeight','Bold')
shg