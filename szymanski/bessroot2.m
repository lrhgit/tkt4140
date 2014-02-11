function [zj0, zj1,niter] = bessroot2(s,rtol)
% Secant-method
b0 = (s - 0.25)*pi;
z00 = b0 + 0.125/b0;
z01 = z00*1.001;
j00 = besselj(0,z00);
test = 1.0;
sit1 = 0;
while (test > rtol)
    sit1 = sit1 + 1;
    j01 = besselj(0,z01);
    dz = - j01*(z01 - z00)/(j01 - j00);
    z = z01 + dz;
    test = abs(dz/z);
    z00 = z01;
    z01 = z;
    j00 = j01;
end
zj0 = z;
b1 = (s + 0.25)*pi;
z10 = b1 - 0.375/b1;
z11 = z10*1.001;
j10 = besselj(1,z10);
test = 1.0;
sit2 = 0;
while (test > rtol)
    sit2 = sit2 + 1;
    j11 = besselj(1,z11);
    dz = - j11*(z11 - z10)/(j11 - j10);
    z = z11 + dz;
    test = abs(dz/z);
    z10 = z11;
    z11 = z;
    j10 = j11;
end
zj1 = z;
niter = sit1 + sit2;