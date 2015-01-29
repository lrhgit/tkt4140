 function tvalue = kanalyt(r,t)
%
% Analytisk løsning av temperaturfordeling T(r,t) i et
% en kule med radius b der 0 <= r <= b.
% Kula med starttemperatur Tk blir sluppet i vann
% med konstant temperatur Tv.
% Denne løsningen forutsetter av H*b/K = 1
% der H er varmeovergangstallet og 
% K er varmeledningstallet.
% alfa er termisk diffusivitet (cm^2/s)
% Tk er starttemperatur i kula (t = 0)
% Tv er konstant temperatur i vannet.
% Meget langsom konvergens for r nær 0 og samtidig
% t nær 0.
%
global alfa b Tk Tv;
fac1 = pi*0.5/b;
tegn = 1; test = 1; epsi = 1.0e-5;
s = 0; m = -1;
f1 = b*1.0e-5; f2 = b*(1-1.0e-5);
if (r >= f1) & (r <= f2)  
 while (test > epsi)
   m = m + 2;  
   bm = m*fac1;  
   bmr = r*bm;
   temp = exp(-alfa*bm^2*t)/m^2;
   ledd = tegn*sin(bmr)*temp;
   s = s + ledd;
   tegn = -tegn;
   test = abs(temp/s);
 end
 tvalue = Tv + 8*b*(Tk - Tv)*s/(pi^2*r);
 return
end
% r = 0
if r < f1
 while (test > epsi)
   m = m + 2;  
   bm = m*fac1;  
   ledd = tegn*exp(-alfa*bm^2*t)/m;
   s = s + ledd;
   tegn = -tegn;
   test = abs(ledd/s);
 end
 tvalue = Tv + 4*(Tk - Tv)*s/pi;
 return
end
% r = b
if r > f2
 while (test > epsi)
   m = m + 2;  
   bm = m*fac1;  
   ledd = exp(-alfa*bm^2*t)/m^2;
   s = s + ledd;
   test = abs(ledd/s);
 end
 tvalue = Tv + 8*(Tk - Tv)*s/pi^2;
end
