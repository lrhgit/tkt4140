%program f76ana
% Beregner den analytiske løsningen
% av det termiske problemet i fig. 7.6
% med bruk av funksjonen txy
clear
yvec = (0 : 0.05 : 0.95)';
n = length(yvec);
T = zeros(n,1);
x = 0.95;
display(x)
n = length(yvec);
for k = 1:n
    y = yvec(k);
    T(k) = txy(x,y);
end
[yvec T]
