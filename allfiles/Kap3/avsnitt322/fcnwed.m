function dy = fcnwed(x,y)
% Brukt i programmet Wedge
dy = zeros(size(y));
if (x > 0.5)
   b2 = 8.0;
else
   b2 = 2.0;
end
if (x > 0)
   dy(1) = b2*y(2)/x;
else
   dy(1) = 2*y(1);
end
dy(2) = y(1); 
