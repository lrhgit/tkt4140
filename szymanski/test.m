% program test
clear
for s = 1:200
    z0 = j0zero2(s);
    z02 = j0test(s);
    feil = abs((z0 - z02)/z0);
    fprintf('%4.0f %12.3e  \n',s,feil);
end

