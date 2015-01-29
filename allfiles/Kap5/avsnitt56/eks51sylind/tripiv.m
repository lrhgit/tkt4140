function [x,fail] = tripiv(a,b,c,d)

 %  Solution of a linear system of algebraic equations with 
 %  a tridiagoal matrix of coefficients using pivoting.  
 %  Equation no. i :                                      
 %                                                              
 %     a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i),          
 %                                         i = 1,2,...neq      
 %                                                              
 %  ======================== Input ============================   
 %                                                                                
 %  a(1:neq) .. vector . Lower diagoneqal. Element a(1)        
 %                                          is not used           
 %  b(1:neq) .. vector . Main diagonal                        
 %  c(1:neq) .. vector . Upper diagoneal. Element c(neq)         
 %                                          is not used           
 %  d(1:neq) .. vector . Right hand side of the system.       
 %                                                               
 %  ========================= Output ==========================   
 %                                                                
 %  x(1:neq) .. vector . The solution vector                  
 %  fail  ..           . fail = 0  -> successful elimination  
 %                       fail = -1 -> singular matrix  
 
 neq = length(b);
 x = zeros(size(b));
 fail = 0;
 
 %===== Reordering =====
 a(1) = b(1);
 b(1) = c(1);
 c(1) = 0;
 
 %===== Elimination =====
 l = 1;
 for k = 1 : neq
    q  = a(k);
    i = k;
    if(l < neq)
       l = l + 1;
    end
    for j = k + 1 : l
       q1 = a(j);
       if(abs(q1) > abs(q))  
          q = q1;
          i = j;
       end  
    end
    if ( q == 0) 
       fail = - 1;
       return  
    end
    if(i ~= k) 
       q  = d(k);
       d(k) = d(i);
       d(i) = q;       
       q  = a(k);
       a(k) = a(i);
       a(i) = q;       
       q  = b(k);
       b(k) = b(i);
       b(i) = q;       
       q  = c(k);
       c(k) = c(i);
       c(i) = q;       
    end
    for i = k + 1 : l       
       q    = a(i)/a(k);       
       d(i) = d(i) - q*d(k);       
       a(i) = b(i) - q*b(k);       
       b(i) = c(i) - q*c(k);       
       c(i) = 0;       
    end    
 end
 
 %===== Backsubstitution =====
 x(neq)   = d(neq)/a(neq); 
 x(neq-1) = (d(neq-1) - b(neq-1)*x(neq))/a(neq-1); 
 for i = neq-2 : -1 : 1    
    q    = d(i) - b(i)*x(i+1);    
    x(i) = (q - c(i)*x(i+2))/a(i);    
 end
 