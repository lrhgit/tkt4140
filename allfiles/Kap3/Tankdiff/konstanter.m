% === konstanter ===
% Beregner konstantene C1,C2,C3 og C4
% som funksjon av beta i den analytiske
% løsningen av vanntank-problemet.
antall = 15;
C = zeros(antall,4);
for k = 1:15
    bt = k;
    J = 4*bt*(cos(bt)^2 + cosh(bt)^2);
    C1 = (2*cos(bt)^2*(bt - 1.0) - sin(2*bt)*bt + bt*(1+exp(-2*bt)))/J;
    C(k,1) = C1;
    C2 = (2*cos(bt)^2*bt + sin(2*bt)*(bt - 1) -(1 + bt)*(1+exp(-2*bt)))/J;
    C(k,2) = C2;
    C3 = (sin(2*bt)*bt + (1 + bt)*2*cos(bt)^2 + bt*(1 + exp(2*bt)))/J;
    C(k,3) = C3;
    C4 = -(2*cos(bt)^2*bt - sin(2*bt)*(1 + bt) +(1 - bt)*(1+exp(2*bt)))/J;
    C(k,4) = C4;
end
beta = [1:15]';
plot(beta,C(:,4));
fprintf(' %5.2f   %14.5e    %14.5e    %14.5e    %14.5e  \n',[beta C]');

% wa = exp(bt*x).*(C1*cos(bt*x) + C2*sin(bt*x))...
%    + exp(-bt*x).*(C3*cos(bt*x) + C4*sin(bt*x)) -(1 - x); 
% [x w wa]
