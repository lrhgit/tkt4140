function thetav2 
% === Ikke-stasjonær varmetransport i konisk kjøleribbe ===
% Som thetav1, men bruker en modifisert versjon av tdma
% Oppgave 1 og 2  i Matlabøving 6, 2011
% Løser også øving 9, 2011, med bruk av theta-skjemaet
% Bruker fletta funksjoner og trenger Matlab ver. 7.x
% theta = 0:   FTCS-skjemaet
% theta = 1/2: Crank-Nicolsons metode
% theta = 1:   Laasonen-skjemaet 
% Ligningen er gitt ved :
%
%      x*T''(x) + 2*T'(x) - b^2*T(x) = x*dT/dt
%          med  0 < x < 1
% Randbetingelser : T(1) = 1 , 2*T'(0) = b^2*T(0)
% der b = beta er Biot-tallet.
% Tidskritt: dt = D*h^2
% Ingen plotting.
% Lokale funksjoner : solut 
%============================================================
clear
N = 9; % Antall deler = N + 1
h = 1/(N + 1); % Skrittlengde
beta = 2; % Biot-tallet
beta2 = beta^2; 
D = 0.4; % Numerisk diffusjonstall
dt = D*h^2;
fac = 1/(1 + beta2*h/2 + (beta2*h)^2/12);
nmax = 60; % Antall tidskritt
tid = nmax*dt;
disp('      *****************************************');
disp('      *        KONISK KJØLERIBBE              *');
disp('      *    theta = 0:   FTCS-skjema           *');
disp('      *    theta = 1:   Laasonen-skjema       *');
disp('      *    theta = 1/2: Crank-Nicolson-skjema *');
disp('      *****************************************');
disp('');
fprintf('\n Antall tidskritt............... %4.0f\n',nmax);
fprintf(' Skrittlengde................... %7.3f\n',h);
fprintf(' Numerisk diffusjonstall........ %7.3f\n',D);
fprintf(' Biot-tall...................... %7.3f\n',beta);
fprintf(' Tidskritt...................... %12.3e\n',dt);
fprintf(' Foreløpt tid................... %12.3e\n\n',tid);
%
% --- Allokerer vektorer 
a = zeros(N,1); b = a; c = a; d = a;
Told = a; Tnew = a; T1 = a; T2 = a; T3 = a;
%
% === theta = 0 : FTCS ===
theta = 0; Told = zeros(N,1);
T1 = solut(theta);
%
% === theta = 1: Laasonen ===
theta = 1; Told = zeros(N,1);
T2 = solut(theta);
%
% === theta = 1/2 : Crank-Nicolson ===
theta = 0.5; Told = zeros(N,1);
T3 = solut(theta);
%
%=== Utskrift av T1, T2 og T3 ===
T1zero = T1(1)*fac; T1 = [T1zero; T1; 1];
T2zero = T2(1)*fac; T2 = [T2zero; T2; 1];
T3zero = T3(1)*fac; T3 = [T3zero; T3; 1];
x = [0:h:1]';
fprintf('      x       T(ftcs)      T(laasonen)   T(crank-n)   \n\n')
fprintf('  %7.3f  %12.4e  %12.4e  %12.4e \n',[x T1 T2 T3]');
%
% === Lokal funksjon solut løser theta-skjemaet ===
%
function vec = solut(theta)

% Da venstre side av ligningsystemet er tids-
% uavhengig , bruker vi en modifisert versjon
% av tdma som vist i lign. 13 i appendiks 9.
for j = 1:N
    a(j) = D*theta*(j - 1);
    b(j) = - (j*(2*D*theta + 1) + D*beta2*h*theta);
    c(j) = D*theta*(j + 1);
end
%=== Eliminasjon,lign. (13a),appendiks 9
c(1) = c(1)/b(1);
for j = 2 : N
    b(j) = b(j)- a(j)*c(j-1);
    c(j) = c(j)/b(j);
end
for n = 1: nmax
    % === Beregning av høyre side 
    for j = 1:N       
        temp1 = 0.0; 
        if j > 1
            temp1 = D*(1-theta)*(1-j)*Told(j-1);
        end
        temp2 = 0;
        if j < N
            temp2 = -D*(1-theta)*(j+1)*Told(j+1);
        end
        d(j) = temp1 +(D*(1-theta)*(2*j + beta2*h)-j)*Told(j)+ temp2;
    end
    % === Innsetting , lign. 13b og c 
    d(N) = -D*(N + 1) + d(N);
    d(1) = d(1)/b(1);
    for j = 2:N
        d(j) = (d(j) - a(j)*d(j-1))/b(j);
    end
    Tnew(N) = d(N);
    for j = N-1:-1 :1
        Tnew(j) = d(j) - c(j)*Tnew(j+1);
    end
    Told = Tnew; % Neste tidskritt
end
vec = Tnew;
end % function solut
end % function combo
