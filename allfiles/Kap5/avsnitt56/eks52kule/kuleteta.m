function kuleteta
% === Eksempel 5.2 med bruk av teta-skjemaet ===
% Beregner temperaturfordeling i et
% en kule med radius b som avkjøles i vann
% Starttemperatur = Tk
% Temperatur i vann  = Tv = konstant
% Beregning ved bruk av FTCS-skjemaet.
% r-skritt dr og det numeriske Fourier-tallet D
% gis og tidskrittet dtau beregnes.
% Stopptiden tauend gis.
%
% Bruker 2.ordens bakoverdifferanser
% i randbetingelsen for r = b = 5cm
% rb = 1 : Bruker randbetingelse med separat ligning
%          for r = 0. Stabil for D <1/3
% rb = 2 : Bruker 2.ordens foroverdifferanser  for r = 0.
%          Stabil for D < 1/2
%     
% Den analytiske løsningen beregnes av
% funksjonen kanalyt. Betingelse for bruk av
% den analytiske løsningen : H*b/K = 1
% der H = varmeovergangstallet og
% K = varmeledningstallet
% 
% Using nested functions (Matlab ver. 7.x)
% teta = 0:   FTCS-scheme
% teta = 1/2: Crank-Nicolson 
% teta = 1:   Laasonen
%
% Output of u. No  plotting.
% Local functions : soltemp, tdma 
% External functions : kanalyta (analytical soltempion)
%============================================================
clear
clear global alfa b Tk Tv;
global alfa b Tk Tv;
b = 5; % radius av kule
m = 50; % Antall deler
alfa = 0.04; % Termisk diffusivitet (cm^2/s)
dr = b/m;     % r-skritt i cm
D = 0.40; % Numerisk Fourier-tall
dtau = D*dr^2/alfa; % Tidsinkrement i sekund
tauend = 600;     % Stopptid i sekunder
ntot = round(tauend/dtau); % antall tidskritt
K = 0.1; % Varmeledningstall (W/cm/C)
H = 0.02 ; %Varmeovergangstall (W/cm^2/C)
dra = dr*H/K;
%Ta = zeros(m+1,1); % Analytiske verdier
Tk = 300; % Starttemperatur i kula
Tv = 20; % Temperatur i vannet.
Told = Tk*ones(m+1,1);% Initialverdier
Tnew = Told;
rb = 2;
N = 50; % No of parts
h = 1/N; % Length of r-step
D = 0.4; % Numerical diffusion number
dt = D*h^2;
nmax = 500; % No. of time-steps
tid = nmax*dt;
disp('      *****************************************');
disp('      *        Varmeledning i kule            *');
disp('      *    teta = 0:   FTCS-skjema            *');
disp('      *    teta = 1:   Laasonen               *');
disp('      *    teta = 1/2: Crank-Nicolson         *');
disp('      *****************************************');
disp('');
fprintf('\n No. of time-steps..............  %6.0f\n',nmax);
fprintf(' r-step length.................. %7.3f\n',h);
fprintf(' Diffusion-number D............. %7.3f\n',D);
fprintf(' Timestep....................... %12.3e\n',dt);
fprintf(' Elapsed time................... %12.3e\n\n',tid);
%
% --- Allocating vectors 
a = zeros(N,1); b = a; c = a; d = a;
Tnew = a;
wa = zeros(N + 1,1);
Told = wa;
r = (0:h:1)';

% === teta = 0 : FTCS ===
teta = 0; Told = 1 - r.^2;
w1 = soltemp(teta);

% === teta = 1: Laasonen ===
% teta = 1; Told = 1 - r.^2;
% w2 = soltemp(teta);
% 
% % === teta = 1/2 : Crank-Nicolson ===
% teta = 0.5; Told = 1 - r.^2;
% w3 = soltemp(teta);

% === Analytical solution ===
% for l = 1: length(r)
%     x = r(l);
%     wa(l) = fcnwa(x,tid);
% end
% === Output of w1, w2 , w3 and wa ===
w1 = [w1; 0]; w2 = [w2; 0]; w3 = [w3; 0];
w1 = 1 - r.^2 - w1; % u-field FTCS
w2 = 1 - r.^2 - w2; % u-field Laasonen
w3 = 1 - r.^2 - w3; % u-field Crank-Nicolson
wa = 1 - r.^2 - wa; % u-field analytical
fprintf('      r       u(ftcs)      u(laasonen)   u(crank-n)    analyt.  \n\n');
s1 = ' %7.3f  %12.4e  %12.4e  %12.4e  %12.4e\n';
for k = 1:length(r)
    fprintf(s1,r(k),w1(k),w2(k),w3(k),wa(k));
end

% === Lokal funksjon soltemp løser teta-skjemaet ===
%
function vec = soltemp(teta)
    tmp1 = D*(1 - teta);   
    for n = 1: nmax
        for j = 2:N
            fac = 1/(j - 1);
            a(j) = D*teta*(1 - fac);
            b(j) = -(1 + 2*D*teta) ;
            c(j) =  D*teta*(1 + fac);            
            d(j) = -tmp1*((1 - fac)*Told(j-1) -tmp1*(1 + fac)*Told(j+1))...
                   - (1 - 2*tmp1)*Told(j);
        end
        b(1) = -(1 + 2*D*teta);
        c(1) =  2*D*teta;
        d(1) = (2*tmp1-1)*Told(1) - 2*tmp1*Told(2);
        Tnew = tdma(a,b,c,d);
        Told = [Tnew;0]; % Updating for the next time-step
    end
    vec = Tnew;
end % function soltemp
    function x = tdma(a,b,c,d)
        m = length(b);
        x = zeros(size(b));
        for k = 2:m
            q = a(k)/b(k-1);
            b(k) = b(k) - c(k-1)*q;
            d(k) = d(k) - d(k-1)*q;
        end
        q = d(m)/b(m);
        x(m) = q;
        for k = m-1 :-1 :1
            q = (d(k) - c(k)*q)/b(k);
            x(k) = q;
        end
    end % function tdma
end % function startup
