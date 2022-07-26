function [g]=G(eta,E,xi,F,d,b1,b2)
a = 2*(1-eta^2)/(pi*E);

X1A = (d^2-b2*xi)/(d*abs(b2-xi));
X1B = (d^2-b1*xi)/(d*abs(b1-xi));

%arrendondamento por erro numerico
if X1A < 1
    X1A =1;
elseif X1B < 1
    X1B =1;
end

bd2 = b2/d;
bd1 = b1/d;

%arrendondamento por erro numerico
if abs(bd2)>1
        bd2=round(bd2);
end
if abs(bd1)>1
        bd1=round(bd1);
end

X1 = (b2-xi)*acosh(X1A)-(b1-xi)*acosh(X1B)+sqrt(d^2-xi^2)*(asin(bd2)-asin(bd1));

xi2 = -1*xi;

X2A = (d^2-b2*xi2)/(d*abs(b2-xi2));
X2B = (d^2-b1*xi2)/(d*abs(b1-xi2));

%arrendondamento por erro numerico
if X2A < 1
    X2A =1;
elseif X2B < 1
    X2B =1;
end

X2 = (b2-xi2)*acosh(X2A)-(b1-xi2)*acosh(X2B)+sqrt(d^2-xi2^2)*(asin(bd2)-asin(bd1));
G1 = a*X1*F;
G2 = a*X2*F;
g = G1 + G2;
end