function [value,isterminal,direction] = Events(t,y)
% Called by program golf45
value = y;
isterminal = [0; 1; 0 ;0];
direction =  [0; -1; 0; 0];