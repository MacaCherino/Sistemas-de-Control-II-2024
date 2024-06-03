% Control sin Observador

clear all; clc;
%% Parámetros del sistema.
a=0.07; be=5; c=150 ; omega=9; 
%Tiempo de muestreo, paso de integración, tiempo de simulación:
Ts=0.1 ; VecesEuler=100; At=Ts/VecesEuler; T=70 ; KMAX=T/Ts; 
%Definición e inicialización de variables
alfa_p=0; fi_p=0; fi_pp=0; h_p=0; ii=1;
t=0:At:T; 
%Matrices del sistema
Mat_Ac=[-a a 0 0; 0 0 1 0; (omega^2) -(omega^2) 0 0 ; c 0 0 0]; %Matriz de  estados.
Mat_Bc=[0; 0; (omega^2)*be; 0]; %Matriz de entrada
Mat_C=[0 0 0 1; % Dos variables de salida, la altura h y fi
 0 1 0 0] 
Mat_D=[0 ; 
 0];
%% Discretizacion del sistema
sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,Mat_D);
sys_d=c2d(sys_c,Ts,'zoh'); 
Mat_A=sys_d.a; 
Mat_B=sys_d.b;
%Matriz de controlabilidad
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B] 
rango=rank(Mat_M)
%Polos de lazo abierto
auto_val=eig(Mat_A);
%% ----------------Cálculo del controlador----------------------------------
Q=diag([100 40 1 0.00012]);
R=10;
[K,P,E]=dlqr(Mat_A,Mat_B,Q,R);
%Ganancia de prealimentación
G=inv(Mat_C(1,:)*inv(eye(4)-Mat_A+Mat_B*K)*Mat_B);
%Polos del controlador a lazo cerrado
aut_controlador=abs(eig(Mat_A-Mat_B*K))

%% Valor de Referencia
%h_ref=-100; % Referencia de altura
h_ref=-100;
%Condiciones iniciales;
h_inic=-500; % Altura inicial
%h_inic=500;
alfa(1)=0; fi(1)=0; fi_p(1)=0; h(1)=h_inic; u(1)=0; ref(1)=h_ref; 
%Vector de estados
x=[alfa(1); fi(1); fi_p(1); h(1)];
x_hat=[0;0;0;0];
%% Simulación
for ki=1:KMAX
    ref(ki)=h_ref; 
    % Ley de control
    u1(ki)=-K*x+G*ref(ki); 
    % Integracion Euler
    for kii=1:Ts/At
        u(ii)=u1(ki);
        ref(ii)=ref(ki);
        alfa_p = a*(fi(ii)-alfa(ii));
        fi_pp = -(omega^2)*(fi(ii)-alfa(ii)-(be*u(ii)));
        h_p = c*alfa(ii);
        alfa(ii+1) = alfa(ii)+At*alfa_p;
        fi_p(ii+1) = fi_p(ii)+At*fi_pp;
        fi(ii+1) = fi(ii)+At*fi_p(ii);
        h(ii+1) = h(ii)+At*h_p;
        ii=ii+1;
    end
    x=[alfa(ii-1); fi(ii-1); fi_p(ii-1); h(ii-1)];
    y_sal=Mat_C*x; 
end
%% Plots
u(ii)=u1(ki);
ref(ii)=ref(ki);

if h_ref >0
    if h_inic > 0
        color = [1 0 0];
    else
        color = [0 0 1];
    end
else
    if h_inic > 0
        color = [1 0 1];
    else
        color = [0.4660 0.6740 0.1880];
    end
end

fz=15;
figure(1);

subplot(2,2,1);
plot(t,alfa,'Color', color,'Linewidth',1.2);
grid on;
title(' $\alpha$ , Angulo de vuelo [rad]', 'Interpreter','latex','FontSize', fz);

subplot(2,2,2);
plot(t,fi,'Color', color,'Linewidth',1.2);
grid on;
title('$\phi$ , Angulo de cabeceo [rad]', 'Interpreter','latex','FontSize', fz);

subplot(2,2,3);
plot(t,fi_p,'Color', color,'Linewidth',1.2);
grid on;
title('$\phi_p$ , Velocidad de cabeceo [rad/s]', 'Interpreter','latex','FontSize', fz);

subplot(2,2,4);
plot(t,h,'Color', color,'Linewidth',1.2);
grid on; hold on;
plot(t,ref,'--k','Linewidth',1);
title('$h$ , Altitud del avion [m]', 'Interpreter','latex','FontSize', fz);

set(gcf,'Color', 'w');

figure(2);
plot(t,u,'Color', color,'Linewidth',1.2);grid on; title('$u$ , Accion de control', 'Interpreter','latex','FontSize', fz); 
set(gcf,'Color', 'w');

