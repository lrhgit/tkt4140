% program storb
% Avsnitt 2.4.1
% Beregner analytisk løsning av en sirkulær 
% vanntank der veggtykkelsen er konstant.
% Bruker her den tilnærmede løsningen
% for store beta-verdier
clear
R = 8.5; % Tankradius
H = 7.95; % Høyde
t = 0.35; % Veggtykkelse 
ny = 0.2;  % Poissons tall

b = (3*(1 - ny^2)*H^4/(R*t)^2)^0.25;
fprintf('beta = %12.5e \n',b);
% s2b =sin(2*b); cb = cos(b); cb2 = cb^2;

% --- Beregner og skriver ut w, w' w'' og w''' ---
%     wp og dwp er bidraget fra partikulærløsningen

dx = 0.1; kmax = 1 + round(1/dx);
w = zeros(kmax,1); dw = w; d2w = w; d3w = w;

for k = 1:kmax 
    x = (k - 1)*dx;
    wp = -(1 - x);
    sbx = sin(b*x); cbx = cos(b*x);
    wh = exp(-b*x)*(cbx + (1 - 1/b)*sbx);
    w(k) = wh + wp;
    dw(k)= 1 + exp(-b*x)*((1 - 2*b)*sbx - cbx );
    d2w(k) = 2*b*exp(-b*x)*(b*sbx - (b-1)*cbx);
    d3w(k) = 2*b^2*exp(-b*x)*((2*b - 1)*cbx - sbx) ;
end
x = (0 : dx :1.0)';
% s1 = '    x         w            w''(x)';
% s2 = '            w''''(x)         w''''''(x)\n \n';
s1 = '    x         w            w''(x)';
s2 = '        m(x)=-w''''(x)   v(x)=-w''''''(x)\n \n';
fprintf([s1,s2]);
fprintf('%6.3f  %13.5e  %13.5e   %13.5e  %13.5e \n',[x w dw -d2w -d3w]');
