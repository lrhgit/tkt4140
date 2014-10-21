%      program unirot

%
%     r  = f'(0) , s = f"(0)
%
   global odetol u2u1 etamax ksimax

%    ===== I n p u t  =====

%    write (nout,9990)
   fprintf('\n ==== I n p u t ==== \n\n')
   u2u1 =   0.5;
   etamax = 4.0;
   ksimax = 5.0;
   rstart = 0.3;
   rend =   0.9;
   rstep =  0.05;
   sstart = 0.1;
   send =   0.4;
   sstep =  0.05;
   fprintf( '   Velocity ratio u2/u1  = %10.3e \n',u2u1);
   fprintf( '   Max. value of eta     = %10.3e \n',etamax);
   fprintf( '   Max. value of ksi     = %10.3e \n',ksimax);
   fprintf( '   Starting value of r   = %10.3e \n',rstart);
   fprintf( '   Ending value of r     = %10.3e \n',rend);
   fprintf( '   Value of r-step       = %10.3e \n',rstep);
   fprintf( '   Starting value of s   = %10.3e \n',sstart);
   fprintf( '   Ending value of s     = %10.3e \n',send);
   fprintf( '   Value of s-step       = %10.3e \n',sstep);

   odetol = 1.0e-5;
   fztol =  1.0e-5;
   fprintf('\n   Max. error in ode45 = %10.3e \n',odetol);
   fprintf('   Max. error in fzero = %10.3e \n',fztol);

   nr = round((rend - rstart)/rstep);
   rend =  rstart + nr*rstep;
   ns = round((send -sstart)/sstep);
   send  = sstart + ns*sstep;
   options = optimset('TolX',fztol);
   sphi = []; rphi = []; spsi = []; rpsi = [];
   fprintf('\n ===  O u t p u t === \n');
   irphi = 0; irpsi = 0;
   for i = 0:ns
       s = sstep*i + sstart;
       rold = rstart;
       phiold = fcnphi(rold,s);
       psiold = fcnpsi(rold,s);
       for j = 1:nr
           r =  rstep*j + rstart;
           phinew = fcnphi(r,s);
           psinew = fcnpsi(r,s);
           if(phinew*phiold <= 0 ) 
%                fprintf('phi: %12.5e  %12.5e %12.5e \n',s,rold,phiold);
%                fprintf('phi: %12.5e  %12.5e %12.5e \n',s,r,phinew);
               bexact = fzero(@fcnphi,[rold r],options,s);
               irphi = irphi + 1;
               sphi(irphi) = s;
               rphi(irphi) = bexact;
           end
           if(psinew*psiold <= 0)
%                fprintf('psi: %12.5e  %12.5e %12.5e \n',s,rold,psiold);
%                fprintf('psi: %12.5e  %12.5e %12.5e \n',s,r,psinew);
               bexact = fzero(@fcnpsi,[rold r],options,s);
               irpsi = irpsi + 1;
               spsi(irpsi) = s;
               rpsi(irpsi) = bexact;
           end
           phiold = phinew ;
           psiold = psinew;
           rold  = r;
       end
   end
   hold on
   if (irphi ~= 0) 
       plot(rphi,sphi);
       fprintf('Zeros in phi : \n'); 
       for i = 1:irphi
           fprintf( ' %13.6e %13.6e \n',rphi(i),sphi(i));
       end
   else
       fprintf('No zeroes found for function phi !\n');
   end
   if (irpsi ~= 0)
       plot(rpsi,spsi);
       fprintf('Zeros in psi : \n'); 
       for i = 1:irpsi
           fprintf( ' %13.6e %13.6e \n',rpsi(i),spsi(i));
       end
   else
       fprintf('No zeroes found for function psi !\n');
   end
   hold off
   grid 

%      1       11x,'* * * * * * * * * * * * * * * * * * * * * * * * *'/
%      2       11x,'*                                               *'/
%      3       11x,'*          Uniform mixing of layers             *'/
%      4       11x,'*                                               *'/
%      5       11x,'*      Zero-search with ode45 and fzero         *'/
%      6       11x,'*                                               *'/
%      7       11x,'* * * * * * * * * * * * * * * * * * * * * * * * *'/)
% 
