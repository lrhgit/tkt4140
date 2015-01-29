% program test
% Motstrømsvarmeveksler
% Tilfelle 1 i avsnitt 6.7 og appendiks 11 
% Plotter i 2d
clear; close
clearvars -global b1 b4 c4 ;
global b1 b4 c4;
N = 100; % antall xskritt
nmax = 60; % antall tidskritt
% Vektoren printvec bestemmer for hvilket tidstrinn
% det blir laget plott.F.eks. dersom vi ønsker
% plott for tidskritt 5,10,20, nmax, settes
% printvec = [5 10 20 nmax]
printvec = [10 25 45];
nstep = 0; % Antall skritt for økning av u0
       % nstep = 0 setter direkte u0 = 1        
if nstep >= nmax
    disp('-- nmax må være > enn nstep ! --');
    return
end
% stepcase = 1 : Lineær økning av u0 
% stepcase = 2 : Økning etter 2.gradspolynom
% stepcase = 3 : Økning etter 3.gradspolynom
% Verdien av stepcase vilkårlig for nstep = 0
tic
stepcase = 1;
dx = 1/N;
u = zeros(N,1); d1 = u; d2 = u;
v = zeros(N,1);
Tyr = 100; Tir = 30;
Tdiff = Tyr - Tir;
alfy = 0.9;
alfi = 0.2;
% b = wi/wy;
% tilfelle 1 : b = 1
b = 1.0; % 
b1 = 1 + 4*N/alfy;
b5 = 4*N/alfy - 1;
b4 = -(1 + 2*N*(b + 1)/alfi);
c6 = 1 - 2*N*(b + 1)/alfi;
c4 = 2*N*(b - 1)/alfi - 1;
b6 = 2*N*(b - 1)/alfi + 1;
dt = dx; % For størst nøyaktighet
t = 0;
x = (0:dx:1)';
% valg = 1;
% hold on;
% ylim([20 100]);
for n = 1: nmax
    t = t + dt;
    if n <= nstep
        f1 = n/nstep;
        f2 = (n+1)/nstep;
        % === Lineær økning av u0 ===
        if stepcase == 1
            u0n   = f1;
            u0np1 = f2;
        end
        % === Økning etter 2.gradspolynom ===
        if stepcase == 2
            u0n   = f1*(2 - f1);
            u0np1 = f2*(2 - f2);
        end
        % === Økning etter 3.gradspolynom ===
        if stepcase == 3
            u0n   = f1*(3 - f1^2)/2;
            u0np1 = f2*(3 - f2^2)/2; 
        end
    else
        u0n = 1;
        u0np1 = 1;
    end
    d1(1) = -u0np1 + b5*u0n - u(1) + v(2) + v(1);
    d2(1) = - (u0n + u0np1) - u(1) + b6*v(1) + c6*v(2);
    for i = 2:N - 1
        d1(i) = b5*u(i-1) - u(i) + v(i+1) + v(i);
        d2(i) = -u(i) - u(i-1) + b6*v(i) + c6*v(i+1);
    end
     d1(N) = b5*u(N - 1) - u(N) + v(N);
     d2(N) = -u(N) - u(N - 1) + b6*v(N);
    [u,v] = bitris2(d1,d2);
%     if (n == printvec(valg))
%         uvec =[1;u];
%         %vvec = [v;0];
%        Ty = Tdiff*uvec + Tir;
%        %Ti = Tdiff*v + Tir;
%        fprintf('Tid = %8.4f \n\n',t);
%        plot(x,Ty);
%        valg = valg + 1;
%        if valg > length(printvec)
%            valg = valg -1;
%        end
%     end  
end
toc
% grid
% title('Temperatur T{_y}','FontSize',14);
% xlabel('x','FontSize',14);
% ylabel('T','FontSize',14,'Rotation',0);
% hold off;
% u =[1;u];
% v = [v;0];
% Ty = Tdiff*u + Tir;
% %Ti = Tdiff*v + Tir;
% x = (0:dx:1)';
% plot(x,Ty);
% grid
%     


