% program CDfig
clear
x = logspace(0.1,8,200)';
n = length(x);
y2 = zeros(n,1);
cd1 = zeros(n,1);
for k = 1 : n
    Re = x(k);
    y2(k) = CDkule2(Re);
end
loglog(x,y2,'k','LineWidth', 1);
grid
% set(gca,'Xtick',[1 10 1.0e2 1.0e3 1.0e4 1.0e5 1.0e6 ]);
shg
