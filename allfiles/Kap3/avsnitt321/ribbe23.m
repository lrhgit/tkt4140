% Program ribbe23
% Beregner eksemplet i avsnitt 3.2.1
% , tilfelle II, med nummerering som i fig.3.2a
% Bruker tdma.
clear all; close all; clc;
set(0,'DefaultLineLineWidth',2,'DefaultAxesFontName','Arial','DefaultAxesFontSize',20);

beta = 5;
n=5;  % Antall ukjente på innsiden av %beregningsområdet (dvs. uten verdier på x(0)på venste rand %og x(N+1) på høyre rand
h=1/(n+1);

a =  ones(n,1); % Underdiagonal

fac = (beta*h)^2;
b = -(2 + fac)*a;    % Hovediagonal
c = a;                        % Overdiagonal

%% Sett inn for grensebetingesler ved x(h)
b(1) = -(2 + 3*fac); 
c(1) = 2.0;

d = zeros(n,1); % Hoeyrre side
d(n) = -1;    % <= theta(n+1)=1   

%% Loes lineart lign system
thet = tdma(a,b,c,d);

%% compute the value at x0
thet0 = (4*thet(1) - thet(2))/3;
thet = [thet0;thet; 1];   %vektor med verdier for x(0), x(h)...x(n+1) dvs med verdiene på randen => n+2 elt


%% analytisk lsning
xa=linspace(0,1,n+2)'; %% ta med x-verdier paa randen => gir n+2 elt 
%xa = (0 : h: 1.0)';
theta = cosh(beta*xa)/cosh(beta);
fprintf('  x      theta     analyt \n');
fprintf('%4.2f  %9.5f  %9.5f \n',[xa thet theta]');

%% Plot resultat
hl=plot(xa,thet,xa,theta,'-.');
h3=legend('v1','a');