function [Sop] = newman_sop_const(s_f,Smax,alpha,R)
A0 = (0.825 - 0.34*alpha + 0.05*alpha^2)*(cos(pi*Smax/(2*s_f)))^(1/alpha);
A1 = (0.415 - 0.071*alpha)*Smax/s_f;
A3 = 2*A0 + A1 - 1;
A2 = 1 - A0 - A1 - A3;
if R >= 0
CF = A0+A1*R+A2*R^2+A3*R^3; %CF = Sop/Smax
Sop = Smax*CF;
else
CF = A0+A1*R;
Sop = Smax*CF;
end
end