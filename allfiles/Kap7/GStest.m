% Program GStest
Ggs = [0 1/2 0 0; 0 1/6 1/3 0; 0 1/18 1/9 1/3 ; 0 1/36 1/18 1/6];
c = [ -1/2; 7/6 ; 49/18; 49/36 ];
x = zeros(4,1);
for k = 1:12
    x = Ggs*x + c;
    fprintf(' %6.4f %6.4f  %6.4f  %6.4f \n',x')
end
    