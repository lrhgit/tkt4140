function [ans,flag,hused] = kuttam(a,b,eta,f,hstart,hmin,epsr)
    % kuttam computes the value y(b) of the differential-equation
    % y'(x) = f(x,y), y(a) = eta using the Kutta-Merson's method.
    % hstart : the initial value of the steplenght h
    % hmin   : the minimum allowed value of h
    % epsr   : relative accuracy for each iteration step.
    % flag   : = 1 if integration is successful;
    %          = -1 if the required accuracy is not attained.
    % ans    : = y(b) if flag = 1, else unchanged
    % hused  : value of the least h used
    eps32 = epsr/32;
    flag = -1;
    x = a; y = eta; h = hstart; hused = hstart;
    last = false;
    % --- verfify that we do not step past x = b
    % 1 continue
    if (x + h >= b)
        last = true;
    end
    if last
        h = b - x;
    end
    k1 = f(x,y);
    % 2 continue
    k2 = f(x + h/3,y + h*k1/3);
    k3 = f(x + h/3,y + h*(k1 + k2)/6);
    k4 = f(x + h/2,y + h*(k1 + 3*k3)/8);
    k5 = f(x + h, y + h*(k1 - 3*k3 + 4*k4)/2);
    d = h*abs(2*k1 - 9*k3 + 8*k4 -k5)/30;
    ytemp = y + h*(k1 + 4*k4 + k5)/6;
    % --- verify that the required accuracy is attained
    if ( d > epsr*abs(ytemp)) % go to 3
        h = h/2;
        if (h < hmin)
            return
        end
    end
    x = x + h;
    y = ytemp;
    % --- finished if last = true
    if last % go to 4
        ans = y;
        flag = 1;
        return
    end
    % test for possible increase of h
    if ( d <= eps32*abs(y))
        h = 2*h;
        % goto 1
    end
    % --- accuracy not attained: half the step
    % 3 continue
    h = h/2;
    if( h < hmin) % go to 5
        return
    end
    if (h < hused)
        hused = h;
    end
    % --- h is h/2: at least two steps left
    last = false;
    % goto 2
    % --- finished
    % 4 continue
    ans = y;
    iflag = 1;
    % 5 continue
    
    
    
    
    