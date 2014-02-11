%program ex69ana
% Beregner den analytiske løsningen
% av eksempel i avsnitt 6.9.
% u(x,t) = cos(t)*sin(x)
% der 0 <= x <= pi
clear; close
x = linspace(0,pi,51);
M = length(x);
u0 = sin(x);
%tell = 1;
for t = 0 : pi/8 : pi/2;
    u = cos(t)*u0;
    %subplot(7,1,tell)
    plot(x,u)
    hold on
    %axis([-5 5 -1 3])
    %tell = tell + 1;
end
grid

        
