% program tcond
% Avsnitt 2.4.1 ---Sylindrisk tank med konstant veggtykkelse---
% Undersøker tilstandstallet ved numerisk beregning av 
% konstantene C1-C4 som funksjon av beta.
clear
% R = 8.5; % Tankradius
% H = 7.95; % Høyde
% t = 0.35; % Veggtykkelse 
% ny = 0.20;  % Poissons tall
% 
% b = (3*(1 - ny^2)*H^4/(R*t)^2)^0.25;
b = 20
% --- Analytiske verdier ---
s2b =sin(2*b); sb = sin(b); cb = cos(b); cb2 = cb^2;
exm2 = exp(-2*b);exp2 =exp(2*b);

J = 4*b*(cb2 + cosh(b)^2);
C1a = (2*cb2*(b-1) - s2b*b + b*(1+exm2))/J;
C2a = (2*cb2*b + s2b*(b-1) - (1+b)*(1+exm2))/J;
C3a = (2*cb2*(1+ b) + s2b*b + b*(1+exp2))/J;
C4a = -(2*cb2*b -s2b*(1+ b) + (1-b)*(1+exp2))/J;
fprintf(' %12.4e  %12.4e   %12.4e  %12.4e \n',C1a,C2a,C3a,C4a);
%
% --- numerisk beregning ---
a = zeros(4,4); c = zeros(4,1); d = c;
a(1,1) = 1; a(1,2) = 0; a(1,3) = 1; a(1,4) = 0;
a(2,1) = 1; a(2,2) = 1; a(2,3) = -1; a(2,4) = 1;
a(3,1) = -sb; a(3,2) = cb; a(3,3) = exm2*sb; a(3,4) = -exm2*cb;
a(4,1) = -(sb + cb); a(4,2) = cb - sb;
a(4,3) = exm2*(cb - sb); a(4,4) = exm2*(sb + cb);
d(1) = 1; d(2) = -1/b; d(3) = 0; d(4) = 0;
c = a\d;
cond(a)
% x = (0 : dx :1.0)';
% s1 = '    x         w            w''(x)';
% s2 = '            w''''(x)         w''''''(x)\n \n';
% fprintf([s1,s2]);
% fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw d2w d3w]');
