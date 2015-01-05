% === tankex ===
%
% Løser lign. i eksempel 2.6 i kompendiet
% med bruk av differanser. 
% Bruker Richardson ekstrapolering.
% Ligning :
%  w''''(x) + 4*b^4*w(x) = -4*b^4*(1 - x), 0 < x <1
% w(0) = 0, w'(0) = 0
% w''(1) = 0, w'''(1) = 0
%

clear
h1 = 0.005; h2 = h1/2; % Skrittlengder
h = [ h1; h2]; % skrittlengde-vektor
N = 1.0./h; % antall intervall og ligninger
W = zeros(N(2),2);
wf = zeros(N(1),1);
% if mod(1,h) > eps
%     error(' N ikke heltall ! ');
% end
bt = 5.0; % parameter beta

for n = 1:2
    neq = N(n);
    ht = h(n);
    fac = 4*(bt*ht)^4;

    % --- Initialisering ---
    e = zeros(neq,1); f = e; a = e; b = a; c = a; d = a;
    w = a; 

    % --- Første ligning ---
    b(1) = 7.0 + fac;
    c(1) = - 4.0;
    f(1) = 1.0;
    d(1) = - fac*(1.0 - ht);

    for k = 2 : neq - 1
        e(k) = 1.0;
        a(k) = - 4.0;
        b(k) = 6.0 + fac;
        c(k) = - 4.0;
        f(k) = 1.0;
        d(k) = - fac*(1 - k*ht);
    end
    b(neq-1) = 5.0 + fac;
    c(neq-1) = -2.0;
    e(neq) = 2.0;
    a(neq) = -3.0;
    b(neq) = fac;
    d(neq) = 0;
    w = penta(e,a,b,c,f,d);
    W(1:neq,n) = w;
end
a = W(:,1); b = W(:,2);
%a = [0;a]; b = [0; b]; w = [0;w];
W = [a b]

% === Richardson ektrapolering ===
%
for k = 1: N(1)
    m = 2*k - 1;
    wf(k) = W(m,2) + (W(m,2) - W(k,1))/3;
end
% === Analytisk løsning ===

J = 4*bt*(cos(bt)^2 + cosh(bt)^2);
C1 = (2*cos(bt)^2*(bt - 1.0) - sin(2*bt)*bt + bt*(1+exp(-2*bt)))/J;
C2 = (2*cos(bt)^2*bt + sin(2*bt)*(bt - 1) -(1 + bt)*(1+exp(-2*bt)))/J;
C3 = (sin(2*bt)*bt + (1 + bt)*2*cos(bt)^2 + bt*(1 + exp(2*bt)))/J;
C4 = -(2*cos(bt)^2*bt - sin(2*bt)*(1 + bt) +(1 - bt)*(1+exp(2*bt)))/J;
x = [h:h(1):1]';
wa = exp(bt*x).*(C1*cos(bt*x) + C2*sin(bt*x))...
   + exp(-bt*x).*(C3*cos(bt*x) + C4*sin(bt*x)) -(1 - x); 
rerr = abs((wf - wa)./wa);
x = [0;x]; wf = [0;wf]; wa = [0;wa]; rerr = [0;rerr];
[x wf wa rerr]
