function [value,isterminal,direction] = social_trigger(t,y,p,flags)
value = y(22)-p.Mc;
isterminal = 1;
direction = 0;
end