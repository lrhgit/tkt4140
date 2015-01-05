% Program shearzsec
% === Jeffrey-Hamel strømning ===
% Bruker sekantmetoden som nullpunktsløser
% Beregner den verdien av produktet Re*alfa som gir
% null skjærspenning ved veggen( df(1)/deta = 0)
% ved å løse ligningen 
% G(m) = K(m) -sqrt(alfa*(alfa + Re/3)) = 0
% K(m) er det fullstendige elliptiske integralet 
% av 1. slag. alfa i radianer.
% K(m) beregnes av Matlabfunksjonen ellipke(m)
% m = 1/(1+D) der D = 1+6*alfa/Re
% som gir m = Re/(2*Re + 6*alfa).
% Lager tabell for 5 < alfag < 90 der
% alfag er i grader : alfa = alfag*pi/180
% Re beregnes da ved å løse ligningen ovenfor.
clear
alfag = 20; % vinkel i grader

% --- Sekantmetoden ---
Re0 = 590/alfag; Re1 = Re0*1.05;
alfa = alfag*pi/180; % Konverterer til radianer
m0 = 0.5*Re0/(Re0 + 3*alfa);
G0 = ellipke(m0) - sqrt(alfa*(alfa + Re0/3));
% --- Startverdier for iterasjonen ---
itmax = 10; epsi = 1.0e-5; it = 0; dRe = 1.0;
while (abs(dRe) > epsi) & (it < itmax)
    it = it + 1;
    m1 = 0.5*Re1/(Re1 + 3*alfa);
    G1 = ellipke(m1) - sqrt(alfa*(alfa + Re1/3));
    dRe = - G1*(Re1 - Re0)/(G1 - G0);
    Re = Re1 + dRe;
    Re0 = Re1;
    Re1 = Re;
    G0 = G1;
    fprintf('%10d %12.5f %12.3e \n',it,Re,dRe);
end

