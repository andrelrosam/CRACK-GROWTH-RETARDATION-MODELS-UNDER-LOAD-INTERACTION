%modelo de willenborg generalizado ou Gallergh
clc
clear
%input carregamento
Pmax = 11.75e3; %N
Rol = [1.375 1.685 2.6];
% P0 = Rol*Pmax; %Overload (N)
R = 0.06;

%input geoemetrico
t = 2.35e-3; %m
b = 50e-3; %m
a0 = 10e-3; %m

%input material
s_0 = 338; %tensao de escoamento (MPa) por Lei Xiao
Kic = 24.17; %tenacidade a fratura (MPam^1/2)
m = 2.84; %constante de walker por Rodrigo
C = 6.3e-11; %coeficiente de walker (m/ciclo) por Rodrigo
g_walker = 0.5; %constante do 2024T3 por Rodrigo
n = 1- g_walker; %constante
DeltaK_th = 1.7; %Utilizado por J.A. Moreto

%calculo do carregamento
smax = Pmax/(2*b*t)*10^-6; %MPa
smin = R*smax; %MPa

%plano de deformacao ou de tensao
if t>= 2.5*(Kic/s_0)^2
    beta = 6; %plano de deformaçao
else 
    beta = 2; %plano de tensao
end

for i=1:1:size(Rol,2)
u = 0; %constante auxiliar 
r=1; %auxiliar
P0 = Rol(i)*Pmax; %Overload (N)
s0 = P0/(2*b*t)*10^-6; %MPa
     
%caculo de propagacao de trinca por ciclo
a=a0;
[Kmax,~,DeltaK] = int_tensao(a,b,smax,smin);
da_dN = C*(DeltaK/(1-R)^n)^m; % m/ciclo
N=1; %numero de ciclos
ai=a+da_dN;
ry0 = 1/(beta*pi)*(Kmax/s_0)^2;
delay = 0;%para dar inicio ao loop

while ai<0.02 %criterio de parada, neste caso e o tamanho de trinca
%calculo do fator de intensidade de tensao
if  ai >= 0.015 && u == 0
    [K0,Kmin,DeltaK] = int_tensao(ai,b,s0,smin);
    ry0i = 1/(beta*pi)*(K0/s_0)^2;
    dadn_before = da_dN;
    u = 1;
else 
    [Kmax,Kmin,DeltaK] = int_tensao(ai,b,smax,smin);
    ry0i = 1/(beta*pi)*(Kmax/s_0)^2;
end

%criterio de retardo
if (ai+ry0i)>(a+ry0)
    a=ai;
    ry0 = ry0i;
else
    j = 0;
end

if j == 0
    delay = delay + 1;
    j = 1;
    Da = ai-a;
    Kreq = K0*sqrt(1-Da/ry0) - Kmax;
    alpha = (1-DeltaK_th/DeltaK)/(Rol(i)-1);
    Kr = alpha*Kreq;
else
    Kr = 0;
end

if (Kmin-Kr) < 0
    Reff = 0;
else
    Reff = (Kmin-Kr)/(Kmax-Kr);
end

da_dN = C*(DeltaK/(1-Reff)^n)^m;

if delay == 1
    dadn_after = da_dN;
end

vetor_a(r) = [ai];
vetor_Reff(r) = [Reff];

ai = ai + da_dN;
N = N+1;
r = r + 1;
end

a_0 = b-P0/(2*t*s_0*10^6);
if a_0>ai
    disp('falha frágil')
else
    disp('falha por escoamento, refazer')
end

data_delay(i) = [delay];
data_dadn_bef(i) = [dadn_before];
data_dadn_aft(i) = [dadn_after];
data_Kol(i) = [K0];
data_N(i) = [N];

figure(1)
plot(vetor_a*1000,vetor_Reff,'LineWidth',2)
hold on
title('Gallagher Model')
xlim([14.5 16.5])
ylim([0 0.1])
xlabel('Crack length - c [mm]')
ylabel('Effective Ratio - R_{eff}')
grid on
legend('R_{ol} = 1.375','R_{ol} = 1.685','R_{ol} = 2.6')
end