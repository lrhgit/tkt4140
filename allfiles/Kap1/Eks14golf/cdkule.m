function cd = cdkule(re)
% Computes the drag-coefficient of a sphere as a function
% of the Reynolds number re. Curve fitted after fig. A-56
% in J.B. Evett & Cheng Liu :
% "Fluid Mechanics & Hydraulics", Schaum's Solved Problems,
%  McGraw-Hill 1989.

if (re <=0.0)
   cd = 0.0;
   return
end
if (re > 8.0e6) 
   cd = 0.2;
   return
end
x0 = log(re)/log(10.0);
if (re> 0.0 & re <=0.5)
   cd = 24.0/re;
   return
end
if (re> 0.5 & re <=100.0)
    p = [4.22 -14.05 34.87 0.658];
    cd = polyval(p,1.0/re);
    return
end
if (re > 100.0 & re <=1.0e4)
    p = [-30.41 43.72 -17.08 2.41];
    cd = polyval(p,1.0/x0);
    return
end
if (re > 1.0e4 & re <=3.35e5)
    p = [-0.1584 2.031 -8.472 11.932];
    cd = polyval(p,x0);
    return
end
if (re > 3.35e5 & re <=5.0e5)
    x1 = log(re/4.5e5)/log(10.0);
    cd = 91.08*x1^4 + 0.0764;
    return
end
if (re > 5.0e5 & re <=8.0e6)
    p = [-0.06335 1.1905 -7.332 14.93];
    cd = polyval(p,x0);
    return
end

