%======================== CDkule ==========================
% Computes the drag-coefficient of a sphere as a function
% of the Reynolds nunber Re. Curve fitted after fig. A-56
% in "Fluid Mechanics & Hydraulics", Schaum's Solved Problems
%=============================================================

function CD = CDkule(Re)
a2=[0.65841 34.8671 -14.046 4.22396];
a3=[2.41329 -17.0825 43.7187 -30.414];
a4=[11.9321 -8.47162 2.03087 -0.15843];
a5=[14.9294 -7.33205 1.19054 -6.33556e-2];
if (Re==0.0)
   CD = 0.0;
   return
end
 if (Re> 8.0e6) 
   CD=0.2;
   return
end
 x0 = log(Re)/log(10.0);
if (Re> 0.0 & Re<=0.5)
   CD=24.0/Re;
   return
end
if (Re> 0.5 & Re<=100.0)
    CD=parab(a2,1.0/Re);
    return
end
if (Re> 100.0 & Re<=1.0e4)
    CD=parab(a3,1.0/x0);
    return
end
if (Re> 1.0e4 & Re<=3.35e5)
    CD=parab(a4,x0);
    return
end
if (Re> 3.35e5 & Re <=5.0e5)
    x1 = log(Re/4.5e5)/log(10.0);
    CD=91.082*x1^4 + 0.0764;
    return
end
if (Re> 5.0e5 & Re<=8.0e6)
    CD=parab(a5,x0);
    return
end
 %-----------------------
 function val = parab(y,x) 
 val = y(1)+x*(y(2)+(y(3)+ y(4)*x)*x);

