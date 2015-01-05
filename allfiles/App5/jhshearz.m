function jhshearz
% --- Jeffrey-Hamel strømning----
% Bruker Matlabfunksjonen fzero som nullpunktsløser.
% Beregner den verdien av produktet Re*alfa som gir
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
close
options = optimset('TolX',1.0e-7);
% --- Overskrift tabell ---
fprintf(' alfa(grd)        Re          Re*alfa \n');
k = 0;
for alfag = 5:5:85
    k = k + 1;
    % --- startverdier 
    if (alfag < 45)
        Re0 = 590/alfag;
    elseif (alfag < 65)
        Re0 = -0.274*alfag + 22.1;
    else
        Re0 = -0.17*alfag + 15.3;
    end
    alfa = alfag*pi/180; % Konverterer til radianer
    Re(k) = fzero(@fcngm,Re0,options); 
    Realfa(k) = Re(k)*alfa;
    fprintf('%10.2f  %12.5e   %12.5e \n',alfag,Re(k), Realfa(k));
end
length(Re)
algv = (5:5:85)';
plot(algv,Realfa);
FW = 'FontWeight'; FS = 'FontSize';
xlabel('\alpha\circ \rightarrow',FS,14);
ylabel('Re\cdot\alpha',FS,14);
title('Re\cdot\alpha with zero shear at the wall',FS,14);
grid
function gval = fcngm(Re)
m = 0.5*Re/(Re + 3*alfa);
gval = ellipke(m) - sqrt(alfa*(alfa + Re/3));
end
end % jhshearz
