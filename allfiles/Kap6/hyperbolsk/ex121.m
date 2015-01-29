%program ex121
clear; close
x = -5 : 0.01 : 5;
tell = 1;
for t = 0 : 0.5 : 3;
    x1 = x + t;
    x2 = x - t;
    for k = 1: 1001
       u(k) = 0.5*(fcnex121(x1(k)) + fcnex121(x2(k)));
    end
    subplot(7,1,tell)
    plot(x,u)
    hold on
    axis([-5 5 -1 3])
    tell = tell + 1;
end

        
