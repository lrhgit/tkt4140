function [s0,s1,etainf] = svalue(beta)
% Used by the program fsksec
% For a given value of beta, two estimates of
% f''(0), s0 and s1 are computed, together with an
% estimate of etainf, after checking for a valid value of beta
% The formulas are given in appendix 3, part 3 of the compendium.
%
p1 = [0.633 -1.68 5.76] ; p2 = [36.76  2.0 5.87];
p3 = [0.125  -0.9 5.463]; betasep = -0.19883768;
% === Compute an estimate of s0 and s1 ===
if(beta <= 1.0)
   s0 = (1.27*(0.2 + beta))^0.56;   % betasep <= beta <= 1.0
   s1 = (1.23*(beta - betasep))^0.54;
else
   s0 = 1.23*(beta^0.454);  % 1.0 < beta <= 1.999
   s1 = (-0.0693*beta + 0.661)*beta + 0.642;
end
% === Compute an estimate of etainf, giving a value
% === of f''(0) around 1.0e-5
if ( beta <= 1.0)
   if ( beta >= 0.0 ) % 0 <= beta <= 1.0
      p = polyval(p1,beta);
      etainf = round(10.0*p)*0.1;
   else
      p = polyval(p2,beta); % betasep <= beta < 0
      etainf = round(10.0*p)*0.1;
   end
else
   p = polyval(p3,beta); % 1 < beta <= 1.999
   etainf = round(10.0*p)*0.1;
end




