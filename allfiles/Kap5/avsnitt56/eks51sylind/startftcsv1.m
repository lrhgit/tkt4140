%========================== startftcsv1 ========================
% === Ikke-stasjonær oppstart i rør, med bruk av FTCS-skjemaet.
% (Szymanskis problem)
% 
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
% Plotter hastighetsprofilet u for hvert tidskritt
%============================================================
clear; clf;
N = 10; % Antall deler
h = 1/N; % Skrittlengde 
D = 0.45; % Numerisk diffusjonstall
dt = D*h^2;
nmax = 80; % Antall tidskritt
disp('      ************************************');
disp('      *           Szymanski              *');
disp('      *          FTCS-skjema             *');
disp('      ************************************');
disp('');
fprintf('\n Antall deler................... %4.0f\n',N);
fprintf(' Antall tidskritt............... %4.0f\n',nmax);
fprintf(' Numerisk diffusjonstall........ %7.3f\n',D);
fprintf(' Tidskritt...................... %12.3e\n',dt);

% --- Allokerer vektorer samt startverdier
r = (0:h:1)';
u = zeros(N + 1,1); 
wnew = 1 - r.^2;
us = wnew;
wold = wnew;
xlim([0 1]);
hold on
for n = 1: nmax
    tid = n*dt;
    for j = 2:N
        fac = 0.5/(j - 1);
        a = D*(1 - fac);
        b = 1 - 2*D;
        c = D*(1 + fac);
        wnew(j) = a*wold(j-1) + b*wold(j) + c*wold(j+1);
    end
    wnew(1) = (1 - 4*D)*wold(1) + 4*D*wold(2);
    wold = wnew;
    u = us - wnew;
    plot(u,r)
end
hold off
FS = 'FontSize';
title('Oppstart av strømning i rør',FS,14);
xlabel('u',FS,14);
ylabel('r',FS,14,'Rotation',0);
% fprintf(' Foreløpt tid................... %12.3e\n\n',tid);
% fprintf('      r         w   \n\n')
% fprintf('  %7.3f  %12.4e \n',[r wnew]');
% plot(r,wnew)
shg;
