% Program beuler
% Løser modell-problem for testing av stabilitet
% av baklengs Eulers metode
% Ligning : y'(x) = -a*y(x), y(0) = 1
% med løsning : y(x) = exp(-a*x);
clear
a = 1;
xmax = 10.0;
dx = 1.25;
nmax = round(xmax/dx); % Number of time-steps
% Initializing of av vectors
x = zeros(nmax,1); y = x; 
y(1) = 1;  x(1) = 0;
clf;
% === Backward Euler method ===
for n = 1 : nmax - 1
    y(n+1) = y(n)/(1 +  a*dx);
    x(n+1) = dx*n;
end 
    plot(x,y,'k')
    hold on
for dx = [ 1.5 2.0 2.5 ]
    nmax = round(xmax/dx); % Number of time-steps
    % Initializing of av vectors
    x = zeros(nmax,1); y = x; 
    y(1) = 1;  x(1) = 0;
    % === Backward Euler method ===
    for n = 1 : nmax - 1
        y(n+1) = y(n)/(1 +  a*dx);
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