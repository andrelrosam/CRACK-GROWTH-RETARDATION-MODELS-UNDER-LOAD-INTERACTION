%modelo de newman dugdale with crack closure
clc
clear

%input carregamento
Pmax = 10e3; %N
Rol = [1.3 1.5 1.75];
R = 0.1;

%input geoemetrico
t = 1.3e-3; %m
b = 50e-3; %m
W = 2*b; %m
ci = 11e-3; %m

%input material
s_0 = 369; %tensao de escoamento (MPa) por Rodrigo
s_u = 490; %tensao de ultimate (MPa) por Rodrigo
s_f = (s_0+s_u)/2; %tensao flow
E = 71.6; %Modulo de elasticidade (MPa) por Rodrigo
Kic = 24.17; %tenacidade a fratura (MPam^1/2) por NASGRO
Kc = 165; %Fator de intensidade de tensao (MPam^1/2) por Golden
nu = 0.35; %coeficiente de poisson por D.chen

%calculo do carregamento
Smax = Pmax/(2*b*t)*10^-6; %MPa
Smin = R*Smax; %MPa

%momento de overload
a_ol = 16.5e-3; %m

%plano de deformacao ou de tensao
if t>= 2.5*(Kic/s_0)^2 %plano de deformaçao
    eta = nu; %propridade do material
    alpha = 3; %constrain factor
else %plano de tensao
    eta = 0; %propridade do material
    alpha = 1; %constrain factor 
end

%coeficientes de elber modificado
C1=4.65e-10; %coeficiente de fit vindo de forman por Rodrigo
C2=2.8; %coeficiente de fit vindo de forman por Rodrigo
C3=0.8901; %coef from threeshold alpha 1
C4=0.0168; %coef from threeshold
C5=Kc; %critical fracture toughness para 1.27 mm

% inputs iniciais
c = ci; %tamanho da trinca
N = 0; %numero de ciclos
n = 10; %numero de elementos de barra

%% calculo dos primeiros elementos de barra

for i=1:1:size(Rol,2)
    
crack_length(1) = c; %catalogo de tamanho de trinca
da_dN_regist(1) = 0; %catalogo de tax de propag
data_sop(1) = 0; %catalogo de Sop;
Num_cycles(1) = N; %catalago de num de ciclos

u = 0; %variavel auxiliar
k = 1; %variavel auxiliar

Pover = Rol(i)*Pmax; %Overload (N)
Sover = Pover/(2*b*t)*10^-6; %Overload (MPa)
    
%reboot de ciclo
c = ci; %tamanho da trinca
N = 0; %numero de ciclos
n = 10; %numero de elementos de barra
n_zone = 10; %numero inicial de elementos
Li_max = []; %vetor de deformacao plastica

%set the valores iniciais
Li_max = zeros(1,10);
rho_max = 0;
da_dN_min = 100;
Sop=0;
Sop_save=0;

% flags auxiliares
flag1=0;
flag2=0;
flag3=0;
flag4=0;

while c<0.035 %criteiro de parada, neste caso e o tamanho de trinca

if  c >= a_ol && u == 0 %condicao de spike load
    S = Sover;
    u = 1;
    flag1=1; %flag para iniciar nova discretizacao de zona plastica
    min_elem = 0.01*rho; %salvando tamanho do elemento minimo ate o momento.
    Sop_save = (data_sop(k)+data_sop(k-1))/2;
else
    S = Smax;
end

rho = newman_plastic_zone(c,S,W,alpha,s_f); %calculo da zona plastica
d = rho + c;

%Controle de discretizacao da zona plastica
    if flag1==0 && flag2==0
        n=n-n_zone;
        n_zone_old = n_zone;
        Li_max_old = Li_max;
        wi = rho/2*[0.3, 0.2, 0.15, 0.12, 0.09, 0.06, 0.04, 0.02, 0.01, 0.01]; %vetor de metade da largura para os 10 primeiros elementos de barra na zona plastica
        n_zone = 10;
        for i=1:n
            Li_max(n_zone+i)=Li_max_old(n_zone_old+i);
        end
        n=n+n_zone;
        if length(Li_max)>n
        Li_max(n+1:length(Li_max))=[];
        end
    else
        n=(n-n_zone);
        n_zone_old = n_zone;
        Li_max_old = Li_max;
        n_zone = floor(rho/min_elem);
        array=repmat(min_elem,1,n_zone);
        wi = array/2;
        for i=1:n
            Li_max(n_zone+i)=Li_max_old(n_zone_old+i); %correcao de posicao do vetor de deformacao para os elementos da crack wave 
        end
        n=n+n_zone;
        if length(Li_max)>n
        Li_max(n+1:length(Li_max))=[]; 
        end
    end

xi = d - abs(wi(1));

for i=2:n
    if i<=n_zone
        xi(i) = xi(i-1) - ( abs(wi(i)) + abs(wi(i-1))); %definindo as posicoes dos bar elements na zona plastica
    else
        wi(i) = crack_ext(i-n_zone)/2; %crack extension = 2wn
        xi(i) = xi(i-1) - ( abs(wi(i)) + abs(wi(i-1))); %definindo as posicoes dos bar elements na crack wave
    end
end

%mudanca de nomeclatura
wj = wi;
xj = xi;

%% calculo de deformacao plastica na zona plastica

for i=1:n
    Somatorio = 0;
for j=1:n
    [b1,b2,B1,B2] = newman_Bk(xj(j),wj(j),W,d);
    [F] = newman_geo_factor(W,d,b1,b2,B1,B2);
    g(i,j) = G(eta,E,xi(i),F,d,b1,b2);
    if j<=n_zone
    Somatorio = Somatorio + alpha*s_f*g(i,j);
    end
end
fx(i) = newman_f(eta,E,W,d,xi(i));
if i<=n_zone
Li_max(i) = S*fx(i) - Somatorio;
end
end

%% a seguir o metodo de gauss para corrigir as tensoes de cada elemento e recalcular a deformacao plastica

s_i = zeros(1,n_zone);
Li = Li_max;
[s_j,itr] = gauss_seidel(S,fx,Li,alpha,s_i,s_f,g,n_zone);

for i=1:n_zone
    Somatorio = 0;
for j=1:n_zone
    Somatorio = Somatorio + s_j(j)*g(i,j);
end
Li_max(i) = S*fx(i) - Somatorio;
end

%% criterio de finalizacao para retardo criado

if (Sop_save*1.05>Sop && flag4==1) %margem de erro usada de 5%
    flag3=0;
    flag4=0;
elseif Sop>Sop_save && flag2==1
    flag4=1;
end

%% confirmando maxima deformacao plastica nos elementos de barra

if (rho>rho_max || (cmax+rho_max)<=c) || flag3==0 %zona atual saiu totalmente da zona max OU zona atual é a nova maxima
   
    %salvando os valores
    rho_max = rho;
    cmax = c;
    xi_max = xi(1:n_zone);
    wi_max = wi(1:n_zone);
    Li_max_ret = Li_max(1:n_zone);
    n_zone_max = n_zone;
    
    if flag1==1 %flag para auxiliar na parada de discretizacao variavel
        flag2=1; %permite discretizacao
        flag3=1; %flag para auxiliar parada de retardo
    else
        flag2=0; %sinaliza fim da discretizacao continua
    end
    
else %(c+rho)>(cmax+rho_max) || (c+rho)<=(cmax+rho_max), zona atual esta no processo de saida da zona max ou ainda reside totalmente dentro da zona max
    d_max = cmax+rho_max;
    Li_max = conf_deformacao_mod(d_max,xi_max,wi_max,Li_max_ret,n_zone_max,xi,Li_max,n_zone);
    flag1=0;
end

%% inicia o calculo para aplicacao de Smin

if n>(n_zone+1)

Vi = [];
s_i = zeros(1,n);
Li = Li_max;
[s_j,itr] = gauss_seidel(Smin,fx,Li,alpha,s_i,s_f,g,n_zone);

for i=1:n
    Somatorio = 0;
for j=1:n
    Somatorio = Somatorio + s_j(j)*g(i,j);
end
Vi(i) = Smin*fx(i) - Somatorio;

%conferencia de displaciment com deformacao plastica recomendada por Newman
erro = Vi(i)/Li_max(i)-1; %margem de erro necessaria por erros numericos de calculo
    if abs(erro)>=0.05
       if Vi(i)>Li_max(i) && i>n_zone
           s_j(i) = 0;
      elseif s_j(i)<0 && i<=n_zone && Vi(i)>Li_max(i)
          Li_max(i) = Vi(i); %corrigi a deformacao plastica em elementos da zona plastica
       end
    end
end

%% Calculo do Sop
c_o = c - crack;

[Sop] = newman_Sop_mod(s_j,xj,wj,Smin,W,c_o,n_zone);

if Sop <= Smin
   warning('Crack is fully open.');
   Sop = 0;
end

else
    Sop = 0; %estou fazendo isso apenas para os dois primeiros loops, onde praticamente nao afeta o calculo de da_dn
end
%% calcular a nova taxa de propagaçao, para computar os ciclos

[da_dN] = newman_elber_cg(S,d,c,W,Sop,C1,C2,C3,C4,C5);

Delta_c = 0.05*rho;

crack = 0;
N_old = N;
N = 0;

%contabilidade de ciclos ate a trinca atingir o passo
while crack < Delta_c && N < 300
crack = crack + da_dN;
N = N + 1;
end

if crack > Delta_c %correcao caso tenha passado o passo maximo permitido
    crack = Delta_c;
    N=N-1;
end

%% aplicando criterio de unificacao para elementos na crack wake
i=(n_zone+1);
while n>=(n_zone+2) && i<n
        if (2*(wi(i)+wi(i+1)))<=(c+Delta_c-xi(i+1)) %condicao de unificacao do elemento da crack wave por newman
            Li_max(i) = (Li_max(i)*wi(i)+Li_max(i+1)*wi(i+1))/(wi(i)+wi(i+1));
            wi(i) = wi(i)+wi(i+1);
            crack_ext(i-n_zone) = crack_ext(i-n_zone)+crack_ext(i+1-n_zone);
            wi(i+1)=[];
            Li_max(i+1)=[];
            xi(i+1)=[];
            crack_ext(i+1-n_zone)=[];
            n=length(wi);
            for j=i:n
                xi(j) = xi(j-1) - ( abs(wi(j)) + abs(wi(j-1))); %so concerto a coordenda para manter o loop funcionando, caso contrario nao é necessario
            end
            i=(n_zone+1);
        else
            i=i+1;
        end
end

%% add contabil para o novo ciclo
n = n + 1; %add um elemento para dar inicio ao novo ciclo
N = N + N_old; %salva o numero de ciclos totais
c = c + crack; %salva o comprimento de trinca atual

%% concertando as coordenadas e deformacoes dos elementos antes de comecar o prox ciclo 
for i=n:-1:(n_zone+1) %loop para transladar saves de deformacao e espessura no vetor antes do novo ciclo
if i==(n_zone+1)
    crack_ext(i-n_zone) = crack; %salva a ultima extensao de trinca que teve
    Leq = defor_eq(crack,Li_max(1:n_zone),wi); %deformacao plastica equivalente
    Li_max(i) = Leq; %salvar as deformacoes plasticas da crack wave que serao imutaveis
else
    crack_ext(i-n_zone) = crack_ext(i-(n_zone+1));
    Li_max(i) = Li_max(i-1); %salvar as deformacoes plasticas da crack wave que serao imutaveis
end
end

%% saves para plot e analise
k=k+1;
crack_length(k) = c;
data_sop(k) = Sop;
da_dN_regist(k) = da_dN;
Num_cycles(k) = N;
n_registro(k) = n;
def_registro(k) = Li_max(1);
coord(k) = xi(1);
end

%% plots
figure(1)
plot(crack_length*1000,data_sop)
hold on
title('Modified Newman Model')
xlim([14 20])
% ylim([20 50])
ylabel('Crack Opening stress - S_{op} [MPa]')
xlabel('Crack length - c [mm]')
grid on
legend('Case 01','Case 02','Case 03')

display(N)
display(max(n_registro))

figure(2)
semilogy(crack_length*1000,da_dN_regist)
hold on
hold on
title('Modified Newman Model')
xlim([15 19])
% ylim([20 50])
ylabel('$log(\frac{da}{dN})$','Interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0)
xlabel('Crack length - c [mm]')
grid on
legend('Case 01','Case 02','Case 03')

end