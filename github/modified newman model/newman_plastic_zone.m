function [rho] = newman_plastic_zone(c,Smax,W,alpha,s_f)
x = (pi*c)/W;
rho = c*((1/x)*asin(sin(x)*sec((pi*Smax)/(2*alpha*s_f)))-1);
end