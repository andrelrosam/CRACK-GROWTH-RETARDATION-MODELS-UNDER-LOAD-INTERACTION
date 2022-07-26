function [Sop] = newman_Sop_mod(s_j,xj,wj,Smin,W,c_o,n_zone)
Soma=0;
n=length(s_j);
for j=(n_zone+1):(n-1)
    b2 = xj(j) + wj(j);
    b1 = xj(j) - wj(j);
    B2 = sin(pi*b2/W)/sin(pi*c_o/W);
    B1 = sin(pi*b1/W)/sin(pi*c_o/W);
    if abs(B2)>1
        B2=round(B2);
    end
    if abs(B1)>1
        B1=round(B1);
    end
    Soma_i = 2*s_j(j)/pi*(asin(B2) - asin(B1));
    Soma = Soma + Soma_i;
end
Sop = Smin - Soma;
end