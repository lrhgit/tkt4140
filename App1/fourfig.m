% program fourfig
% Plotter funksjonen gm(x) = sin(mx) , 0 <= x <= 2*pi,
% for forskjellige verdier av m, m >=1
clear; close;
x = 0:0.02:2*pi;
y = zeros(length(x));
xlim([-0.5 6.5]);
hold on
y = 0;
plot (x,y,'k');
for m = [1 6]
    y = sin(m*x);
    plot (x,y,'k');
end
grid
axis off
hold off
