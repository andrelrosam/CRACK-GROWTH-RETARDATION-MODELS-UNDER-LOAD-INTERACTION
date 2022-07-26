function [b1,b2,B1,B2] = newman_Bk(xj,wj,W,d)
b2 = xj + wj;
b1 = xj - wj;
B2 = sin(pi*b2/W)/sin(pi*d/W);
B1 = sin(pi*b1/W)/sin(pi*d/W);
%arrendondamento por erro numerico 
if abs(B2)>1
        B2=round(B2);
end
if abs(B1)>1
        B1=round(B1);
end

end