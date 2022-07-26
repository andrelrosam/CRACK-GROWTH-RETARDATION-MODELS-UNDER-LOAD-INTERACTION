function [da_dN] = newman_elber_cg(Smax,d,c,W,Sop,C1,C2,C3,C4,C5)
F = sqrt(sec(pi*d/W));
Kmax = Smax*sqrt(pi*c)*F;
dKop = C3*(1-C4*Sop/Smax);
dKeff = (Smax-Sop)*sqrt(pi*c)*F;
da_dN = C1*dKeff^C2*(1-(dKop/dKeff)^2)/(1-(Kmax/C5)^2);
if da_dN<0 || Sop>Smax %condicao de crack arrest
    da_dN = 5e-8;
end
end