function gval = fcngm(Re)
% Used by program jhshearz and jhshearz3
global alfa;
m = 0.5*Re/(Re + 3*alfa);
gval = ellipke(m) - sqrt(alfa*(alfa + Re/3));
