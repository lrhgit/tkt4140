% ===================== Program jhfigur ===============================
% Solution of Jeffrey-Hamel-flow, appendix 5,
% using a shooting technique with the secant-method as a root-finder.
% Restricted case with alpha^2 << Re*alpha
% The program generates a plot of f for different values of Re*alpha.

% The Jeffrey-Hamel equation in this case is given by :
%  f'''(eta) + 2*Re*alpha*f(eta)*f'(eta)  = 0
%  f(0) = 1, f'(0) = 0, f(1)=0
% In the  program we put y(1)= f, y(2) = f' and y(3) = f'' 
% 
% s = f''(0) = y(3)(0)
%
% Using function fcnjh2
% =====================================================================
clear; clear global Rea;
global Rea 
Reavec = [ -50.0 -10.0 5.0 10.3128];
s0vec  = [ -0.011 -0.545 -3.820 -6.87519];
etaspan = [0 1]; % The integration-interval
% We have found suitable starting-values s0 and s1
% from the program plotjh.
antall = length(Reavec);
options = odeset('RelTol',1.0e-5);
% Start with alpha = 0 : Parabolic profil
% with f = 1 - eta^2;
eta0 = linspace(0,1,50);
f0 = 1 -eta0.^2; 
clf
plot(eta0,f0,'k')
axis([0 1 0 1]);
grid
FS = 'FontSize';
xlabel('\eta',FS,14,'FontWeight','Bold')
ylabel('f',FS,14)
title('Jeffrey-Hamel flow',FS,14)
hold on
for k = 1:antall
    Rea = Reavec(k);
    s0 = s0vec(k);
    s1 = s0*1.05;
    % ---- Compute phi0 ----
    f0 = [1.0 0.0  s0];
    [eta,f] = ode45(@fcnjh2,etaspan,f0);
    phi0 = f(end,1);
    % Starting values for the iteration
    itmax = 10; epsi = 1.0e-5; it = 0; ds = 1;
    % -------- Start of iteration -------------------
    while(abs(ds) > epsi) & (it < itmax)
        it = it + 1;
        f0 = [1.0 0.0 s1];
        [eta,f] = ode45(@fcnjh2,etaspan,f0,options);
        phi1 = f(end,1);
        ds = -phi1*(s1 - s0)/(phi1 - phi0);
        s = s1 + ds;
        s0 = s1;
        s1 = s;
        phi0 = phi1;
        %fprintf('%10d %12.6f %12.3e\n',it,s,ds);
    end
    % -------- Plotting f as a function of eta
    plot(eta,f(:,1),'k')
end