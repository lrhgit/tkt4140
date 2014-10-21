% ================== Program jhshearz3old ============================
%
% --- Jeffrey-Hamel strømning ved avløsning (df(1)/deta = 0) ---
% Beregner f som funksjon av eta direkte med
% bruk av Jacobi-funksjonen sn for 
% for en gitt verdi av vinkelen alfa.
% sn finnes fra Matlabfunksjonen ellipj.
%
% Bruker Matlabfunksjonen fzero som nullpunktsløser
% for å beregne den verdien av produktet Re*alfa som gir
% null skjærspenning ved veggen( df(1)/deta = 0)
% ved å løse ligningen 
% G(m) = K(m) - sqrt(alfa*(alfa + Re/3)) = 0
% G(m) er programert i funksjonen fcngm.
% K(m) er det fullstendige elliptiske integralet 
% av 1. slag. alfa i radianer.
% K(m) beregnes av Matlabfunksjonen ellipke(m)
% m = 1/(1+ D) der D = 1+6*alfa/Re
% som gir m = Re/(2*Re + 6*alfa).
% Lager tabell for 5 <= alfag <= 85 der
% alfag er i grader : alfa = alfag*pi/180
% Re beregnes da ved å løse ligningen ovenfor.
%
clear; clear global alfa;
global alfa
options = optimset('TolX',1.0e-9);
alfag = 20;
% --- startverdier 
if (alfag < 45)
    Re0 = 590/alfag;
elseif (alfag < 65)
    Re0 = -0.274*alfag + 22.1;
else
    Re0 = -0.17*alfag + 15.3;
end
%
%  Beregner den kritiske verdien av Re
%  for den gitte alfag-verdien  
%
alfa = alfag*pi/180; % Konverterer til radianer
Re = fzero(@fcngm,Re0,options); 
Rea = Re*alfa;
fprintf('=== Critical Jeffrey-Hamel flow === \n\n');
fprintf('Angle in degrees (alfag)........ %7.3f \n',alfag);
fprintf('Reynolds number (Re) ............ %13.7e \n',Re);
fprintf('Product Re*alpha ................ %13.7e \n\n',Rea);
fprintf('    eta         f            f'' \n\n');
for k = 1:19
    eta(k) = 0.05*k;
    b = sqrt(3/(alfa*(Re + 3*alfa)));
    ksi = eta(k)/b;
    D = 1.0 + 6*alfa/Re;
    m = 1/(1 + D);
    sn = ellipj(ksi,m);
    f(k) = 1 - sn^2;
    %   Computing df/deta:
    fac = 2*Re*alfa/3; fk = f(k);
    df(k) = - sqrt((1 - fk)*(fac*fk*(fk + 1) + 4*alfa^2*fk));
end
f = [1 f 0]; eta = [0 eta 1]; df = [0 df 0];
fprintf('%7.3f  %13.7e  %13.7e \n',[eta; f; df]);
plot(eta,f,'k',eta,df,'k-.');
grid
st1 = ('Critical case.  ');
st2 = sprintf('Re\\alpha = %5.3f , \\alpha = %3.1f',Rea,alfag);
st = [st1 st2];
xlabel('\eta','FontSize',14,'FontWeight','Bold');
ylabel('f , f''','FontSize',14);
title(st,'FontWeight','Bold')
legend('f','f''')
shg