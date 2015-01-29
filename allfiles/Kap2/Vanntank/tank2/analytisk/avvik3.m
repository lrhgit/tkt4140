% program avvik3
% Denne versjonen er for t1 = 0 -> a = 1
% Sirkulær vanntank med variabel veggtykkelse.
% m0 = dimensjonsløst moment for x = 0
% m0a = asymptotisk dimensjonsløst moment for x = 0
% v0 = dimensjonsløs skjærkraft for x = 0
% v0a = asymptotisk dimensjonsløs skjærkraft for x = 0
% Beregner forskjellen (i %) mellom m0 og m0a og imellom
% v0 og v0a som funksjon av beta for a = 1
% der a = (t0 - t1)/t0 med t1 = 0;
clear
a = 1;
bpoints = 100;
beta = linspace(1.5,5,bpoints)';
rerrm0 = zeros(bpoints,1); rerrv0 = rerrm0;

for k = 1: bpoints
    b = beta(k);
    ro = (b/a)*sqrt(2);
    z0 = 2*ro;
    A = zeros(4,4); B = zeros(4,1); C = B;
   [berz0,beiz0] = berbei(z0);
   [berdz0,beidz0] = berdbeid(z0);
   [kerz0,keiz0] = kerkei(z0);
   [kerdz0,keidz0] = kerdkeid(z0);
   A(1,1) = berdz0; A(1,2) = beidz0;
   A(1,3) = kerdz0; A(1,4) = keidz0;
   A(2,1) = -(2*berdz0 + z0*beiz0);
   A(2,2) = (-2*beidz0 + z0*berz0); 
   A(2,3) = -(2*kerdz0 + z0*keiz0);
   A(2,4) = (-2*keidz0 + z0*kerz0);
   
   A(3,1) = 0;
   A(3,2) = 0;
   A(3,3) = 0;
   A(3,4) = 0;

   A(4,1) = 0;
   A(4,2) = 0;
   A(4,3) = 0;
   A(4,4) = 0;

   fac = 0;
   
   B(1) = 1;
   B(2) = 0;
   B(3) = 0;
   B(4) = 0;
   C(1) = -A(2,2)/(A(1,2)*A(2,1) -A(1,1)*A(2,2));
   C(2) =  A(2,1)/(A(1,2)*A(2,1) -A(1,1)*A(2,2));
   y = 1 ;
   u1 = C(1)*(4*berdz0 -0.5*z0^2*beidz0 + 2*z0*beiz0);
   u2 = C(2)*(4*beidz0 + 0.5*z0^2*berdz0 - 2*z0*berz0);
   u3 = 0;
   u4 = 0;
   m0 = -(0.5*sqrt(y)*(u1 + u2 + u3 + u4) - 2*fac)*a^2;
   u1 = C(1)*(2*beidz0 - z0*berz0);
   u2 = -C(2)*(2*berdz0 + z0*beiz0);
   u3 = 0;
   u4 = 0;
   v0 = a^3*0.5*ro^2*sqrt(y)*(u1 + u2 + u3 + u4);
   fac2 = 4*(4*b -a);
   m0a = (32*b^2*(b-1) + a*(8*b -a))/fac2;
   v0a =- (32*b^2*(2*b^2 - b*(1 + a) + a) -a^2*(8*b - a))/fac2;
   rerrm0(k) = ((m0 - m0a)/m0)*100;
   rerrv0(k) = ((v0 - v0a)/v0)*100;
   %fprintf('%6.3f  %11.4f  %11.4f   %11.4f  %11.4f \n',b, m0, m0a, v0, v0a);
end
% s1 = '    beta       m0            m0a';
% s2 = '           v0         v0a \n \n';
% fprintf([s1,s2]);
% fprintf('%6.3f  %11.2f  %11.2f   %11.2f  %11.2f \n',b, m0, m0a, v0, v0a);
% plot(beta,m0,beta,m0a,beta,v0,beta,v0a)
FS = 'FontSize'; FW = 'FontWeight';
plot(beta,rerrm0,'k',beta,rerrv0,'k')
grid
title('% avvik i m{_0} og v{_0}.  \alpha = 0',FS,14)
xlabel('\beta',FS,14)
ylabel('% avvik',FS,14)
shg