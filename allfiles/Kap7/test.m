% Program GStest
Ggs = [0 1/2 0 0; 0 1/6 1/3 0; 0 1/18 1/9 1/3 ; 0 1/36 1/18 1/6];
c = [ -1/2; 7/6 ; 49/18; 49/36 ];
x = zeros(4,1);
m = 12;
A = zeros(m,4);
for k = 1:m
    x = Ggs*x + c;
    A(k,:) = x;
    fprintf(' %8.6f %8.6f  %8.6f  %8.6f \n',x')
end
    