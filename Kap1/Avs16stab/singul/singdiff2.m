% program singdiff2
% Beregner løsning av singulær diff.ligning
% i avsnitt 1.6.2. 
% Lagrer data på fil for plottingS
clear; 
x = linspace(0,1)';
n = length(x);
A = zeros(n,5);
A(:,1) = x;
epsvec = [0.03 0.1 0.2 0.4];
for k = 2:5
    epsilon = epsvec(k-1);
    fac = exp(1/epsilon) - 1;
    A(:,k) = (exp(x/epsilon)-1)/fac;
end
save figdata.dat A -ascii

% grid on
% xlabel('x','FontSize',14,'FontWeight','Bold')
% ylabel('y','FontSize',14,'FontWeight','Bold','Rotation',0)
% title('Grensesjikt ved x = 1','FontSize',14,'FontWeight','Bold')
% shg