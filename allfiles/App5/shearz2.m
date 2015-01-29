% Program shearz2
%
% === Jeffrey-Hamel strømning ved avløsning ===
% Beregner eta som funksjon av f direkte fra
% det elliptiske integralet for df(1)/deta = 0
% for en gitt verdi av vinkelen alfa.
% Integralet beregnes av funksjonen iellip1
%
% Bruker Matlabfunksjonen fzero som nullpunktsløser
% for å beregne den verdien av produktet Re*alfa som gir
% null skjærspenning ved veggen( df(1)/deta = 0)
% ved å løse ligningen 
% G(m) = K(m) - sqrt(alfa*(alfa + Re/3)) = 0
% K(m) er det fullstendige elliptiske integralet 
% av 1. slag. alfa i radianer.
% K(m) beregnes av Matlabfunksjonen ellipke(m)
% m = 1/(1+ D) der D = 1+6*alfa/Re
% som gir m = Re/(2*Re + 6*alfa).
% Lager tabell for 5 <= alfag <= 85 der
% alfag er i grader : alfa = alfag*pi/180
% Re beregnes da ved å løse ligningen ovenfor.
%
clear
global alfa
options = optimset('TolX',1.0e-9);
% --- Overskrift tabell ---
%fprintf(' alfa(grd)        Re          Re*alfa \n');
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
fprintf('alfag = %7.3f  Re = %13.7e \n\n',alfag,Re);
fprintf('    f          f''           eta              eta/b\n');
for k = 1:9
    f = 0.1*k;
    x = sqrt((1 - f)/f);
    D = 1.0 + 6*alfa/Re;
    kc  = sqrt(D/(1 + D));
    b = sqrt(3/(alfa*(Re + 3*alfa)));
    etab = iellip1(x,kc);
    eta = etab*b;
    %   Beregner df/deta:
    fac = 2*Re*alfa/3;
    df = sqrt((1 - f)*(fac*f*(f + 1) + 4*alfa^2*f));
    fprintf('%7.3f  %13.7e  %13.7e %13.7e\n',f,df,eta,etab);
end