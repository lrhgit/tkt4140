% program tanalyt2
% Avsnitt 2.4.1
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen er konstant.
% Leser inn beta-verdien direkte
clear
beta = input('beta = ? ');
fprintf('     beta = %7.3f\n',beta);
b = beta;
s2b =sin(2*b); cb = cos(b); cb2 = cb^2;

J = 4*b*(cb2 + cosh(b)^2);
C1 = (2*cb2*(b-1) - s2b*b + b*(1+exp(-2*b)))/J;
C2 = (2*cb2*b + s2b*(b-1) - (1+b)*(1+exp(-2*b)))/J;
C3 = (2*cb2*(1+ b) + s2b*b + b*(1+exp(2*b)))/J;
C4 = -(2*cb2*b -s2b*(1+ b) + (1-b)*(1+exp(2*b)))/J;

% --- Beregner og skriver ut w, w' w'' og w''' ---
%     wp og dwp er bidraget fra partikulærløsningen

dx = 0.1; kmax = 1 + round(1/dx);
w = zeros(kmax,1); dw = w; d2w = w; d3w = w;

for k = 1:kmax 
    x = (k - 1)*dx;
    wp = -(1 - x);
    sbx = sin(b*x); cbx = cos(b*x);
    wh = exp(b*x)*(C1*cbx + C2*sbx) + exp(-b*x)*(C3*cbx + C4*sbx);
    w(k) = wh + wp;
    t1 = exp(b*x)*((cbx - sbx)*C1 + (sbx + cbx)*C2);
    t2 = exp(-b*x)*(-(cbx + sbx)*C3 + (cbx - sbx)*C4);
    dwp = 1;
    dw(k) = b*(t1 + t2) + dwp;
    t1 = exp(b*x)*(-C1*sbx + C2*cbx);
    t2 = exp(-b*x)*(C3*sbx - C4*cbx);
    d2w(k) = 2*b^2*(t1 + t2);
    t1 = exp(b*x)*(-C1*(sbx + cbx) + C2*(cbx - sbx));
    t2 = exp(-b*x)*(C3*(cbx - sbx) + C4*(cbx + sbx));
    d3w(k) = 2*b^3*(t1 + t2);
end
x = (0 : dx :1.0)';
s1 = '    x         w            w''(x)';
s2 = '            w''''(x)         w''''''(x)\n \n';
fprintf([s1,s2]);
fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw d2w d3w]');
