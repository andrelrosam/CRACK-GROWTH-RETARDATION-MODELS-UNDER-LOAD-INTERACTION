function [X]=newman_f(eta,E,W,d,xi)
X = 2*(1-eta^2)/E*sqrt((d^2-xi^2)*sec(pi*d/W));
end