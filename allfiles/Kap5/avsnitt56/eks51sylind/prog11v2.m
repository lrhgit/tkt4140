%====================== prog11v2 ====================
% === Ikke-stasjonær oppstart i rør. 
% (Szymanskis problem)
% Øving11, 2007, med bruk av FTCS-skjemaet.
% Lik prog11v1, men kjører til stasjonær løsning.
% Hastighetsfelt : u(r,t) = us(r) -w(r,t)
% Stasjonærløsning : us = 1 - r^2
% Ligningen for w(r,t) er gitt ved :
%
%     dw/dt = w''(r) + w'(r)/r
%          med  0 < r < 1 , w = w(r,t)
%
% Randbetingelse : w(1,t) = w(-1,t) = 0 
% Symmetribetingelse: w'(0,t) = 0
% Startbetingelse : w(r,0) = 1- r^2;
% Stasjonærløsningen for w er 0
% Tidskritt: dt = D*h^2
%=====================================================
clear
N = 2 % Antall deler
h = 1/N; % Skrittlengde 
D = 0.3414; % Numerisk diffusjonstall
dt = D*h^2;
disp('      ************************************');
disp('      *           Szymanski              *');
disp('      *          FTCS-skjema             *');
disp('      ************************************');
disp('');
fprintf('\n Antall deler.................. %4.0f\n',N);
fprintf(' Numerisk diffusjonstall........ %7.4f\n',D);
fprintf(' Tidskritt..................... %12.3e\n',dt);

% --- Allokerer vektorer samt startverdier
r = (0:h:1)'; 
wnew = 1 - r.^2;
wold = wnew;
wtest = 1.0e-4; n = 0; wsum = 1; 
while (wsum > wtest) 
    n = n + 1;
    for j = 2:N
        fac = 0.5/(j - 1);
        a = D*(1 - fac);
        b = 1 - 2*D;
        c = D*(1 + fac);
        wnew(j) = a*wold(j-1) + b*wold(j) + c*wold(j+1);
    end
    wnew(1) = (1 - 4*D)*wold(1) + 4*D*wold(2);
    wold = wnew;
    wsum = abs(sum(wnew)/(N + 1));
    %wmaks = abs(max(wnew));
    if abs(wnew(1)) < abs(wnew(2))
        fprintf('=== Instabilitet ved skritt %7.0f \n',n);
        break
    end
end
tid = n*dt;
fprintf(' Antall tidskritt............... %4.0f\n',n);
fprintf(' Foreløpt tid................... %12.3e\n\n',tid);
fprintf('      r         w   \n\n')
fprintf('  %7.3f  %12.4e \n',[r wnew]');
plot(r,wnew)
shg;
