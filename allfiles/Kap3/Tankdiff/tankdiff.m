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

% --- Initialisering ---
e = zeros(N,1); f = e; a = e; b = a; c = a; d = a;
w = a; 
fac = 4*(bt*h)^4;

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
x = (0:h:1)';
% === Analytisk løsning ===
% J = 4*bt*(cos(bt)^2 + cosh(bt)^2);
% C1 = (2*cos(bt)^2*(bt - 1.0) - sin(2*bt)*bt + bt*(1+exp(-2*bt)))/J;
% C2 = (2*cos(bt)^2*bt + sin(2*bt)*(bt - 1) -(1 + bt)*(1+exp(-2*bt)))/J;
% C3 = (sin(2*bt)*bt + (1 + bt)*2*cos(bt)^2 + bt*(1 + exp(2*bt)))/J;
% C4 = -(2*cos(bt)^2*bt - sin(2*bt)*(1 + bt) +(1 - bt)*(1+exp(2*bt)))/J;
% wa = exp(bt*x).*(C1*cos(bt*x) + C2*sin(bt*x))...
%    + exp(-bt*x).*(C3*cos(bt*x) + C4*sin(bt*x)) -(1 - x); 
%[x w wa];

% === Beregner helningen w'(x) ===
dw = zeros(N+1,1);
dw(1) = 0.0;
for k = 2:N
    dw(k) = 0.5*(w(k+1) - w(k-1))/h;
end;
dw(N+1) = (w(N+1)-w(N))/h;

% === Beregner momentet m(x) ===
m = zeros(N+1,1); 
fac2 = 1/h^2;
m(1) = -2*w(2)*fac2;
m(N+1) = 0;
for k = 2:N
    m(k) = -fac2*(w(k+1) - 2*w(k) + w(k-1));
end;

% === Beregner skjærkrafta v(x). Bruker trapesmetoden ===  
v = zeros(N+1,1);
v(N+1) = 0;
s = 0;
fac3 = -4*bt^4;
for k = N+1:-1:2
    z = 1-(k-2)*h;
    s = s + 0.5*h*(w(k) + w(k-1));
    v(k-1) = fac3*(z^2*0.5 + s);
end;
s1 = '    x         w            w''(x)';
s2 = '        m(x)=-w''''(x)   v(x)=-w''''''(x)\n \n';
fprintf([s1,s2]);
fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw m v]');