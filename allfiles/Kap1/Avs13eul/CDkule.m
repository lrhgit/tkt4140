function CD = CDkule(Re)
% Computes the drag-coefficient of a sphere as a function
% of the Reynolds nunber Re. Curve fitted after fig. A-56
% in Evett & Liu:%  "Fluid Mechanics & Hydraulics",
% Schaum's Solved Problems McGraw-Hill 1989

if (Re <=0.0)
   CD = 0.0;
   return
end
if (Re > 8.0e6) 
   CD = 0.2;
   return
end
x0 = log10(Re);
if (Re> 0.0 & Re <=0.5)
   CD = 24.0/Re;
   return
end
if (Re> 0.5 & Re <=100.0)
    p = [4.22 -14.05 34.87 0.658];
    CD = polyval(p,1.0/Re);
    return
end
if (Re > 100.0 & Re <=1.0e4)
    p = [-30.41 43.72 -17.08 2.41];
    CD = polyval(p,1.0/x0);
    return
end
if (Re > 1.0e4 & Re <=3.35e5)
    p = [-0.1584 2.031 -8.472 11.932];
    CD = polyval(p,x0);
    return
end
if (Re > 3.35e5 & Re <=5.0e5)
    x1 = log10(Re/4.5e5);
    CD = 91.08*x1^4 + 0.0764;
    return
end
if (Re > 5.0e5 & Re <=8.0e6)
    p = [-0.06335 1.1905 -7.332 14.93];
    CD = polyval(p,x0);
    return
end

