function [Kmax,Kmin,DeltaK] = int_tensao(a,b,smax,smin)
%F = (1-0.5*a/b+0.326*(a/b)^2)/sqrt(1-a/b); %fator geometrico
F = sqrt(sec(pi*a/(2*b))); %fator geometrico
Kmax = F*smax*sqrt(pi*a);
Kmin = F*smin*sqrt(pi*a);
DeltaK = Kmax - Kmin;   
end