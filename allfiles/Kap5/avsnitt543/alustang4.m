function alustang4 
% === Ikke-stasjonær varmetransport i aluminiumstang ===
% Tester theta-skjema for sprang i startbetingelser.
% Som alustang3, men bruk av subplot.
%
% theta = 0:   FTCS-skjemaet
% theta = 1/2: Crank-Nicolsons metode
% theta = 1:   Laasonen-skjemaet 
%
% Ligningen er gitt ved :
%     dT/dt = a*T''(x) , 0 < x < 1
%
% a = 100 mm/s^2 (Termisk diffusivitet)
% Stavlengde L = 300 mm
% Tidskritt: dt = 0.25s
% x-skritt = h
% N = L/h - 1; % Antall deler = N + 1
% x(j) = h*j, j = 0,1, ...,N + 1 , 
%
% Startbetingelse for t = 0:
% T = T0 for 0 <= x < 100
% T = Ts for 100 <= x <= 200
% T = T0 for 200 < x <= 300
% 
% Randbetingelser : T(0) = T0 , T(N+1) = T0
% 
% Lokal funksjon : solut
%============================================================
close;
L = 300; % L = 300 mmm
h = 2.5; % dx 
N = L/h - 1; % Antall deler = N + 1
dt = 0.25; % Tidskritt i sekund
alfa = 100; % Termisk diffusivitet i mm^2/s
D = alfa*dt/h^2; % Numerisk Fourier-tall
nmax = 16; % Antall tidskritt
tid = nmax*dt;
theta = 1.0;
%
disp('      *****************************************');
disp('      *          Aluminiumstang               *');  
disp('      *    theta = 1:   Laasonen-skjema       *');
disp('      *    theta = 1/2: Crank-Nicolson-skjema *');
disp('      *****************************************');
disp('');
fprintf('\n Antall tidskritt............... %4.0f\n',nmax);
fprintf(' Skrittlengde................... %7.3f\n',h);
fprintf(' theta-verdi.................... %7.3f\n',theta);
fprintf(' Numerisk diffusjonstall........ %7.3f\n',D);
fprintf(' Diffusivitet................... %7.3f\n',alfa);
fprintf(' Tidskritt...................... %12.3e\n',dt);
fprintf(' Foreløpt tid................... %12.3e\n\n',tid);
%
% --- Allokerer vektorer 
a = zeros(N,1); b = a; c = a; d = a;
Told = d; Tnew = d; Tval = d; Tstart = d;
%
%  --- Startverdier ---
T0 = 20; % grader C
Ts = 270; % grader C 
js = L/(3*h); % L/3 = 100 mm
for k = 1: N
    Tstart(k) = T0;
end
for k = js : 2*js
    Tstart(k) = Ts;
end
x = [0:h:300]';
Told = Tstart;
Tstart = [T0; Tstart; T0];

hold on
plot(x,Tstart)
for n = 1: nmax
    Tval = solut(theta);
    Tval = [T0; Tval; T0];
    plot(x,Tval)
end
grid
hold off
%
% === Lokal funksjon solut løser theta-skjemaet ===
%
function vec = solut(theta)
% Da venstre side av ligningsystemet er tids-
% uavhengig, bruker vi en modifisert versjon
% av tdma som vist i lign. 13 i appendiks 9.
J = 1:N; % Indeksvektor
a(J) = D*theta;
b(J) = - (1 + 2*D*theta); 
c(J) = D*theta;  
%
%=== Eliminasjon,lign. (13a),appendiks 9
c(1) = c(1)/b(1);
for j = 2 : N
    b(j) = b(j)- a(j)*c(j-1);
    c(j) = c(j)/b(j);
end
% === Beregning av høyre side 
temp1 = 1 - 2*D*(1 - theta);
d(1) = - D*T0 -(D*(1 - theta)*Told(2) + temp1*Told(1));
d(N) = - D*T0 -(D*(1 - theta)*Told(N-1) + temp1*Told(N));
for j = 2:N-1
    d(j) = -(D*(1 - theta)*(Told(j-1) + Told(j+1)) + temp1*Told(j));
end
% === Innsetting , lign. 13b og c 
d(1) = d(1)/b(1);
for j = 2:N
    d(j) = (d(j) - a(j)*d(j-1))/b(j);
end
Tnew(N) = d(N);
for j = N-1:-1 :1
    Tnew(j) = d(j) - c(j)*Tnew(j+1);
end
Told = Tnew; % Neste tidskritt
vec = Tnew;
end % function solut
end % function alustang4

