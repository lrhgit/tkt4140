% program avvik
% Avsnitt 2.4.1
% Sirkulær vanntank med konstant veggtykkelse.
% m0 = dimensjonsløst moment for x = 0
% m0a = asymptotisk dimensjonsløst moment for x = 0
% v0 = dimensjonsløs skjærkraft for x = 0
% v0a = asymptotisk dimensjonsløs skjærkraft for x = 0
% Beregner forskjellen (i %) mellom m0 og m0a og imellom
% v0 og v0a som funksjon av beta.

clear
bpoints = 40;
beta = linspace(1.5,5,bpoints)';

rerrm0 = zeros(bpoints,1); rerrv0 = rerrm0;
for k = 1 : bpoints
    b = beta(k);
    s2b =sin(2*b); cb = cos(b); cb2 = cb^2;
    J = 4*b*(cb2 + cosh(b)^2);
    m0 = -4*b^2*(2*cb2 - 2*b*cosh(b)^2 + sinh(2*b) - s2b)/J;
    v0 = -8*b^3*(b*(s2b + sinh(2*b)) + cb2 - cosh(b)^2)/J;
    m0a = 2*b*(b - 1);
    v0a = -2*b^2*(2*b - 1);
    rerrm0(k) = ((m0 - m0a)/m0)*100;
    rerrv0(k) = abs((v0 - v0a)/v0)*100;
end
% s1 = '    beta       m0            m0a';
% s2 = '           v0         v0a \n \n';
% fprintf([s1,s2]);
% fprintf('%6.3f  %11.2f  %11.2f   %11.2f  %11.2f \n',[beta m0 m0a v0 v0a]');
% plot(beta,m0,beta,m0a,beta,v0,beta,v0a)
FS = 'FontSize'; FW = 'FontWeight';
plot(beta,rerrm0,'k',beta,rerrv0,'k')
grid
title('Relativt avvik i m{_0} og v{_0}',FS,14)
xlabel('\beta',FS,14)
ylabel('% avvik',FS,14)
shg