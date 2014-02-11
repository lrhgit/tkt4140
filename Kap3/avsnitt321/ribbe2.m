%============ Program ribbe2 =================
% Lser varmeledningsproblemet i avsnitt 3.2
% der diff.ligningen er gitt i ligning 3.2.4
% Dette er en kombinert versjon av ribbe21 og
% ribbe22 med tdma som lser.
%============================================
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

n=10;  % Antall ukjente på innsiden av 
          % beregningsområdet (dvs. med verdien på x(1)på venste rand  men uten og x(n+1) på høyre rand
h=1/n;
%h = 0.25; % skrittlengde
%n = round(1.0/h); % Antall ligninger

beta = 5.0; % Biottallet er beta^2

fact = -(2.0 + (beta*h)^2);
a = ones(n,1) ; % underdiagonal
b = fact*ones(n,1); % hoveddiagonal
c = a ; % overdiagonal 
d = zeros(n,1) ; % hyre side av ligningsystemet
d(n) = -1.0;

%------------------------------------------
% Tilfelle 2, versjon 1 : lign 2.11 + 2.13
%-----------------------------------------
c(1) = 2.0 ; % for versjon 1
theta = tdma(a,b,c,d); % Ls ligningsystemet
theta(n+1)= 1.0;

%x = (0.0:h:1.0)';
x=linspace(0,1.0,n+1)';

antheta = cosh(beta*x)/cosh(beta); % Analytisk lsning
rerr = abs((antheta - theta)./antheta); % Relativ feil
fprintf(' Tilfelle 2, versjon 1 \n\n')
fprintf('    %6.3f  %10.5f %12.3e \n',[x theta rerr]');

%------------------------------------
% Tilfelle 2 , versjon 2 : lign 2.17
%-----------------------------------
b(1) = 2.0;
c(1) = -( 2.0 - (beta*h)^2);
theta2 = tdma(a,b,c,d); % Ls ligningsystemet
theta2(n+1)= 1.0;
rerr = abs((antheta - theta)./antheta); % Relativ feil
x = (0.0:h:1.0)';
fprintf('\n Tilfelle 2, versjon2 \n\n')
fprintf('    %6.3f  %10.5f %12.3e \n',[x theta rerr]');


%% Plot resultat
hp=plot(x,theta,x,antheta,x,theta2,':');
h3=legend('v1','a','v2');