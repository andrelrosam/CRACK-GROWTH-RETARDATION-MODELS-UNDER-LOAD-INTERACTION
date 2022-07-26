function [F] = newman_geo_factor(W,d,b1,b2,B1,B2)
y = d/W;
bd2 = b2/d;
bd1 = b1/d;
%arrendondamento por erro numerico
if abs(bd2)>1
        bd2=round(bd2);
end
if abs(bd1)>1
        bd1=round(bd1);
end
x = (asin(B2) - asin(B1))/(asin(bd2) - asin(bd1));
F = x*sqrt(sec(pi*y));
end