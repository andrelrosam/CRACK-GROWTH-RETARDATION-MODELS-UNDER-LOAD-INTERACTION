%modelo de wheeler
clc
clear

%input carregamento
Pmax = 11.75e3; %N
Rol = [1.375 1.685 2.6];
R = 0.06;

%input geoemetrico
t = 2.35e-3; %m
b = 50e-3; %m
a0 = 10e-3; %m

%input material
s_0 = 338; %tensao de escoamento (MPa) por Lei Xiao
Kic = 24.17; %tenacidade a fratura (MPam^1/2) por NASGRO
m = 2.84; %constante de walker por Rodrigo
C0 = 6.3e-11; %coeficiente de walker (m/ciclo) por Rodrigo
g_walker = 0.5; %constante do 2024T3 por Rodrigo
n = 1- g_walker; %constante por Rodrigo
C = C0/((1-R)^(m*n)); %coeficeinte de walker (m/ciclo)/(MPam^1/2)
gamma = 3.71; %parametro de fixacao de wheeler para 2024 por Rodrigo

%calculo do carregamento
smax = Pmax/(2*b*t)*10^-6; %MPa
smin = R*smax; %MPa

%plano de deformacao ou de tensao
if t>= 2.5*(Kic/s_0)^2
    beta = 6; %plano de deformacao
else 
    beta = 2; %plano de tensao
end

for i=1:1:size(Rol,2)
u = 0; %variavel auxiliar
p = 1; %variavel auxiliar
data_phir = [];
vetor_a = [];
P0 = Rol(i)*Pmax; %Overload (N)
s0 = P0/(2*b*t)*10^-6; %Overload (MPa)
    
%primeiro ciclo
a=a0;
[Kmax,~,DeltaK] = int_tensao(a,b,smax,smin);
da_dN = C*DeltaK^m; % m/ciclo
N=1; %numero de ciclos
ai=a+da_dN; %incremento de trinca
ry0 = 1/(beta*pi)*(Kmax/s_0)^2; %raio de zona plastica

%para dar inicio ao loop
Kmax_i = 1;
delay = 0;

%calculo de propagacao de trinca por ciclo

while ai<0.02 %criterio de parada, neste caso e o tamanho de trinca
%calculo do fator de intensidade de tensao
if  ai >= 0.015 && u == 0 %condicao de spike load
    [K0,~,DeltaK] = int_tensao(ai,b,s0,smin);
    ry0i = 1/(beta*pi)*(K0/s_0)^2;
    dadn_before = da_dN;
    u = 1;
else %condicao normal
    [Kmax_i,~,DeltaK] = int_tensao(ai,b,smax,smin);
    ry0i = 1/(beta*pi)*(Kmax_i/s_0)^2;
end

%criterio de yield zone
if (ai+ry0i)>(a+ry0)
    a=ai;
    ry0 = ry0i;
else
    j = 0;
end

if j == 0
    delay = delay + 1;
    j = 1;
end

%implementacao final do coeficiente de wheeler
Da = ai-a;
phir = (ry0i/(ry0+a-ai))^gamma;
da_dN = C*DeltaK^m;
da_dNR = phir*da_dN;

data_phir(p) = [phir];
vetor_a(p) = [ai];
data_ciclo(p) = [N];
p = p+1;

%registro de dados abaixo
if delay == 1
    dadn_after = da_dNR;
end

ai = ai + da_dNR;
N = N+1;
end

%verificacao de falha
a_0 = b-P0/(2*t*s_0*10^6);
if a_0>ai
    disp('falha frágil')
else
    disp('falha por escoamento, refazer')
end

figure(1)
plot(vetor_a*1000,data_phir,'LineWidth',2)
hold on
title('Wheeler Model')
xlim([14 16.5])
%ylim([0 0.1])
xlabel('Crack length - c [mm]')
ylabel('Retardation parameter - C_{p}')
grid on
legend('R_{ol} = 1.375','R_{ol} = 1.685','R_{ol} = 2.6')

%salvamento de dados
data_delay(i) = [delay];
data_dadn_bef(i) = [dadn_before];
data_dadn_aft(i) = [dadn_after];
data_Kol(i) = [K0];
data_N(i) = [N];
end