function cd = cdkule2(re)
% Computes the drag-coefficient cd of a sphere as a function
% of the Reynolds nunber re. Fitted Formulae from
% Clift, Grace && Weber : "Bubbles, Drops and Particles",
% Academic Press 1978

if (re <=0.0)
   cd = 0.0;
   return
end
if (re > 1.0e6) 
   cd = 0.19 - 8.0e4/re;
   return
end
x0 = log10(re);
if (re> 0.0 && re <= 0.01)
   cd = 24.0/re + 9/2;
   return
end
if (re> 0.01 && re <= 20.0)
    p = 0.1315*re^(0.82 - 0.05*x0);
    cd =  24*(1 + p)/re;
    return
end
if (re> 20.0 && re <= 260.0)
    cd =  24*(1 + 0.1935*re^0.6305)/re;
    return
end
if (re > 260.0 && re <= 1.5e3)
    p = x0*(-1.1242 + x0*0.1558);
    cd = 10.0^(1.6435 + p);
    return
end
if (re > 1.5e3 && re <= 1.2e4)
    p = x0*(2.5558 + x0*(-0.9295 + x0*0.1049));
    cd = 10.0^(-2.4571 + p);
    return
end
if (re > 1.2e4 && re <= 4.4e4)
    p = x0*(0.6370 - x0*0.0636);
    cd = 10.0^(-1.9181 + p);
    return
end
if (re > 4.4e4 && re <= 3.38e5)
    p = x0*(1.5809 - x0*0.1546);
    cd = 10.0^(-4.3390 + p);
    return
end
if (re > 3.38e5 && re <= 4.0e5)
    cd = 29.78 - 5.3*x0;
    return
end
if (re > 4.0e5 && re <= 1.0e6)
    cd =  0.1*x0 - 0.49;
    return
end

