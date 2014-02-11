% program test2
clear 
antall =2000;
tic
for s = 1:antall
    z0 = j0zero3(s);
end
toc
tic
for s = 1: antall
    z02 = j0zero(s);
end
toc

