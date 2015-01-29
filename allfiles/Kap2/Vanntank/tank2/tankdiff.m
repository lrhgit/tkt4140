% === tankdiff ===
%
% Løser lign. i eksempel 2.6 i kompendiet
% med bruk av differanser.
% Ligning :
%  w''''(x) + 4*b^4*w(x) = -4*b^4*(1 - x), 0 < x <1
% w(0) = 0, w'(0) = 0
% w''(1) = 0, w'''(1) = 0
%
clear

h = 0.01; % skrittlengde
N = 1/h; % antall intervall og ligninger
if mod(1,h) > eps
    error(' N ikke heltall ! ');
end
bt = 5.0; % parameter beta
fac = 4*(bt*h)^4;

% --- Initialisering ---
e = zeros(N,1); f = e; a = e; b = a; c = a; d = a;
w = a; 

% --- Første ligning ---
b(1) = 7.0 + fac;
c(1) = - 4.0;
f(1) = 1.0;
d(1) = - fac*(1.0 - h);

for k = 2 : N - 1
    e(k) = 1.0;
    a(k) = - 4.0;
    b(k) = 6.0 + fac;
    c(k) = - 4.0;
    f(k) = 1.0;
    d(k) = - fac*(1 - k*h);
end
b(N-1) = 5.0 + fac;
c(N-1) = -2.0;
e(N) = 2.0;
a(N) = -3.0;
b(N) = fac;
d(N) = 0;

w = penta(e,a,b,c,f,d);
w = [0; w];
x = [0:h:1]';
% === Analytisk løsning ===
J = 4*bt*(cos(bt)^2 + cosh(bt)^2);
C1 = (2*cos(bt)^2*(bt - 1.0) - sin(2*bt)*bt + bt*(1+exp(-2*bt)))/J;
C2 = (2*cos(bt)^2*bt + sin(2*bt)*(bt - 1) -(1 + bt)*(1+exp(-2*bt)))/J;
C3 = (sin(2*bt)*bt + (1 + bt)*2*cos(bt)^2 + bt*(1 + exp(2*bt)))/J;
C4 = -(2*cos(bt)^2*bt - sin(2*bt)*(1 + bt) +(1 - bt)*(1+exp(2*bt)))/J;
wa = exp(bt*x).*(C1*cos(bt*x) + C2*sin(bt*x))...
   + exp(-bt*x).*(C3*cos(bt*x) + C4*sin(bt*x)) -(1 - x); 
[x w wa]
% theta = [theta; 1];
% % === Beregner theta'(x) fra lign. 2.47) ===
% x = [0:h:1]';
% n = length(theta);
% dtheta = zeros(n,1); % Initialisering
% dtheta(1) = bt2a*theta(1);
% s = 0;
% for k = 2: n
%     s = s + 0.5*h*(theta(k) + theta(k-1));
%     if k > m
%         dtheta(k) = bt2b*s/x(k);
%     else
%         dtheta(k) = bt2a*s/x(k);
%     end
% end
% fprintf('   x       theta       theta'' \n\n');
% x2 = x(m:n);
% theta2 = theta(m:n);
% dtheta2 = dtheta(m:n);
% dtheta2(1) = ka*dtheta(m)/kb;
