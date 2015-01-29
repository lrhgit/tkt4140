function start2 
% === Oppstart av strømning i rør (Szymanskis problem) ===
% Som prog7v1, men bruker modifisert versjon av tdma
% Løser øving 11, 2007, med bruk av theta-skjemaet
% Bruker fletta funksjoner og trenger Matlab ver. 7.x
% theta = 0:   FTCS-skjemaet
% theta = 1/2: Crank-Nicolsons metode
% theta = 1:   Laasonen-skjemaet 
% Hastighetsfelt : u(r,t) = us(r) - w(r,t)
%                  der us = 1 - r^2
% Ligningen for w er gitt ved :
%
%     dw/dt = w'' + w'/r , w = w(r,t)
%          med  0 < r < 1
%
% Randbetingelser : w(1) = w(-1) = 0, w'(0) = 0 
% Startbetingelse : w(r,0) = 1- r^2 = us
% Tidskritt: dt = D*h^2
% Ingen plotting.
% Lokale funksjoner : solut 
%============================================================
clear
N = 50; % Antall deler
h = 1/N; % Skrittlengde
D = 0.4; % Numerisk diffusjonstall
dt = D*h^2;
nmax = 500; % Antall tidskritt
tid = nmax*dt;
disp('      *****************************************');
disp('      *  Impulsiv oppstart av rørstrømning    *');
disp('      *    theta = 0:   FTCS-skjema           *');
disp('      *    theta = 1:   Laasonen-skjema       *');
disp('      *    theta = 1/2: Crank-Nicolson-skjema *');
disp('      *****************************************');
disp('');
fprintf('\n Antall tidskritt............... %4.0f\n',nmax);
fprintf(' Skrittlengde................... %7.3f\n',h);
fprintf(' Numerisk diffusjonstall........ %7.3f\n',D);
fprintf(' Tidskritt...................... %12.3e\n',dt);
fprintf(' Foreløpt tid................... %12.3e\n\n',tid);
%
% --- Allokerer vektorer 
a = zeros(N,1); b = a; c = a; d = a;
wa = zeros(N + 1,1);
wold = a; wnew = wold;
r = (0:h:1-h)';
tic
% === theta = 0 : FTCS ===
theta = 0; wold = 1 - r.^2;
w1 = solut(theta);

% === theta = 1: Laasonen ===
theta = 1; wold = 1 - r.^2;
w2 = solut(theta);

% === theta = 1/2 : Crank-Nicolson ===
theta = 0.5; wold = 1 - r.^2;
w3 = solut(theta);
% === Analytisk løsning ===
r = [r;1];
for l = 1: length(r)
    x = r(l);
    wa(l) = fcnwa(x,tid);
end
toc
%=== Utskrift av u1, u2 , u3 og ua ===
w1 = [w1; 0]; w2 = [w2; 0]; w3 = [w3; 0];
w1 = 1 - r.^2 - w1; % u for FTCS
w2 = 1 - r.^2 - w2; % u for Laasonen
w3 = 1 - r.^2 - w3; % u Crank-Nicolson
wa = 1 - r.^2 - wa; % u analytisk
fprintf('      r       w(ftcs)      w(laasonen)   w(crank-n)    analyt  \n\n')
fprintf('  %7.3f  %12.4e  %12.4e  %12.4e  %12.4e \n',[r w1 w2 w3 wa]');

% === Lokal funksjon solut løser theta-skjemaet ===
%
function vec = solut(theta)

% Da venstre side av ligningsystemet er tids-
% uavhengig , bruker vi en modifisert versjon
% av tdma som vist i lign. 13 i appendiks 9.
tmp1 = D*(1 - theta);
for j = 2:N
    fac = 0.5/(j - 1);
    a(j) = -D*theta*(1 - fac);
    b(j) = (1 + 2*D*theta) ;
    c(j) = - D*theta*(1 + fac);
end
b(1) = 1 + 4*D*theta;
c(1) = - 4*D*theta;
%=== Eliminasjon,lign. (13a),appendiks 9
c(1) = c(1)/b(1);
for j = 2 : N
    b(j) = b(j)- a(j)*c(j-1);
    c(j) = c(j)/b(j);
end
for n = 1: nmax
    % === Beregning av høyre side 
    for j = 2:N       
        fac= 0.5/(j - 1);
        tmp2 = 0;
        if j < N
            tmp2 = (1 + fac)*wold(j+1);
        end
        d(j) = tmp1*((1 - fac)*wold(j-1)+ tmp2) +(1 - 2*tmp1)*wold(j);
    end
    d(1) = wold(1) + 4*tmp1*(wold(2) - wold(1));
    d(N) = tmp1*(1 - fac)*wold(N-1) +(1 - 2*tmp1)*wold(N);
    % === Innsetting , lign. 13b og c 
    d(1) = d(1)/b(1);
    for j = 2:N
        d(j) = (d(j) - a(j)*d(j-1))/b(j);
    end
    wnew(N) = d(N);
    for j = N-1:-1 :1
        wnew(j) = d(j) - c(j)*wnew(j+1);
    end
    wold = wnew; % Neste tidskritt
end
vec = wnew;
end % function solut
end % function start2
