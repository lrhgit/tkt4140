function jhshearz3
%
% --- Jeffrey-Hamel strømning ved avløsning (df(1)/deta = 0) ---
% Beregner f og f' som funksjon av eta direkte med
% bruk av Jacobi-funksjonene sn, cn og dn for 
% for en gitt verdi av vinkelen alfa.
% Matlabfunksjonen ellipj beregner Jacobifunksjonene.
%
% Bruker Matlabfunksjonen fzero som nullpunktsløser
% for å beregne den verdien av produktet Re*alfa som gir
% null skjærspenning ved veggen( df(1)/deta = 0)
% ved å løse ligningen 
% G(m) = K(m) - sqrt(alfa*(alfa + Re/3)) = 0
% G(m) er programmert i funksjonen fcngm.
% K(m) er det fullstendige elliptiske integralet 
% av 1. slag. alfa i radianer.
% K(m) beregnes av Matlabfunksjonen ellipke(m)
% m = 1/(1+ D) der D = 1+6*alfa/Re
% som gir m = Re/(2*Re + 6*alfa).
% Lager tabell for 5 <= alfag <= 85 der
% alfag er i grader : alfa = alfag*pi/180
% Re beregnes da ved å løse ligningen ovenfor.
%
close
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
b = sqrt(3/(alfa*(Re + 3*alfa)));
D = 1.0 + 6*alfa/Re;
m = 1/(1 + D);
eta = zeros(1,19); f = eta; df = eta;
for k = 1:19
    eta(k) = 0.05*k;
    ksi = eta(k)/b;
    [sn,cn,dn] = ellipj(ksi,m);
    f(k) = 1 - sn^2;
    df(k) = -2*sn*cn*dn/b;% Computing df/deta:
end
f = [1 f 0]; eta = [0 eta 1]; df = [0 df 0];
fprintf('%7.3f  %13.7e  %13.7e \n',[eta; f; df]);
plot(eta,f,'k',eta,df,'k-.');
grid
st1 = ('Critical case.  ');
st2 = sprintf('Re\\alpha = %5.3f , \\alpha = %3.1f',Rea,alfag);
st = [st1 st2];
xlabel('\eta','FontSize',14,'FontWeight','Bold');
ylabel('f , f ''','FontSize',14);
title(st,'FontWeight','Bold')
legend('f','f ''')
function gval = fcngm(Re)
m = 0.5*Re/(Re + 3*alfa);
gval = ellipke(m) - sqrt(alfa*(alfa + Re/3));
end % fcngm
end % jhshearz3
