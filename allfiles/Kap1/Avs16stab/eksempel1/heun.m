% Program heun
% Løser modell-problem for testing av stabilitet
% av Heuns metode
% Ligning : y'(x) = -a*y(x), y(0) = 1
% med løsning : y(x) = exp(-a*x);
% Stabilitetsgrense : hmax = 2/abs(a)
clear
a = 1;
xmax = 20.0;
dx = 1.25;
nmax = round(xmax/dx); % Number of time-steps
% Initializing of av vectors
x = zeros(nmax,1); y = x; 
y(1) = 1;  x(1) = 0;
clf;
% === Heun's method ===
for n = 1 : nmax - 1
    y(n+1) = (1 -  a*dx + (a*dx)^2/2)*y(n);
    x(n+1) = dx*n;
end 
    plot(x,y,'k')
    hold on
for dx = [  1.5 1.99  ]
    nmax = round(xmax/dx); % Number of time-steps
    % Initializing of av vectors
    x = zeros(nmax,1); y = x; 
    y(1) = 1;  x(1) = 0;
    % === Heun's method ===
    for n = 1 : nmax - 1
        y(n+1) = (1 -  a*dx+ (a*dx)^2/2)*y(n);
        x(n+1) = dx*n;
    end 
    plot(x,y,'k')
end
grid on
%FS = 'FontSize';
%st = sprintf('\\Deltax = %5.2e',dx);
%title(st,FS,13);
%xlabel('x',FS,14);
%ylabel('y',FS,13);