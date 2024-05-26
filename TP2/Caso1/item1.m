% Control sin Observador
clear all; clc

%% Parámetros del sistema
Laa=4.994e-4; J=2.408e-9;Ra=20;B=0;Ki=9.647e-3;Km=60.529e-3;

%Paso, tiempo de simulación, pasos:
h=5e-7 ; tF=10; pasos=(tF/h);

%Extracción de datos
tabla = readtable("Curvas_Medidas_Motor_2024.xls");
tD = tabla{:,1};
tLD = tabla{:,5};
%Defino el torque
i=0;
tL = zeros(1,pasos);
for t=0:h:0.6
    i=i+1;
    [~,pos] = min(abs(tD-i*h));
    tL(i) = tLD(pos);
    tL(pasos/2+i) = tLD(pos);
end

%% Definición e inicialización de variables
ia_p=0; w_p=0; w_pp=0; theta_p=0; 
t=0:h:tF; ia=0:h:tF; omega=0:h:tF; w_p=0:h:tF; theta=0:h:tF;
u=linspace(0,0,pasos+1); ref=linspace(0,0,pasos+1);
thetaRef=(pi/2); %Referencia de ángulo del motor
Max_U=24; %Tension nominal del motor

tc=5; %Tiempo de cambio de la referencia (5s).
est=0; %Estado de la referencia.

ii=1; kk=0; 
%Condiciones iniciales
ia(1)=0; omega(1)=0; theta(1)=0; wp(1)=0; u(1)=0;
%Matrices del sistema lineal
Mat_A=[-Ra/Laa -Km/Laa 0 ; Ki/J -B/J 0 ; 0 1 0];
Mat_B=[1/Laa; 0 ; 0];
Mat_C=[0 0 1]; % La salida es el ángulo theta.
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B]; %Matríz de controlabilidad
disp('Polos de lazo abierto');
auto_val=eig(Mat_A) %Polos de lazo abierto

%% Cálculo del controlador por asignación de polos
c_ai=poly(Mat_A); %Ecuación característica de lazo cerrado
Mat_W=[c_ai(3) c_ai(2) 1 ; c_ai(2) 1 0 ; 1 0 0]; %Matríz W
Mat_T=Mat_M*Mat_W; %Matríz de transformación
A_controlable=inv(Mat_T)*Mat_A*Mat_T; % Verificación de matriz T

%% Ubicación de los polos de lazo cerrado en mui:
mui(1)=-4e5; mui(2)=-5e4; mui(3)=-6e4; % Sin saturación;
%mui(1)=-2e6; mui(2)=-6e5; mui(3)=-5.3e5; % Con saturación;
alfa_i=poly(mui);
%Controlador
K=(fliplr(alfa_i(2:end)-c_ai(2:end)))*(inv(Mat_T));
disp('Polos de lazo cerrado controlador');
eig(Mat_A-(Mat_B*K)) %Polos de lazo cerrado
Gj=-inv(Mat_C*inv(Mat_A-(Mat_B*K))*Mat_B); %Ganancia de prealimentacion

%%
while(ii<(pasos+1))
 kk=kk+h;
 if(kk>tc) % Cada 5 s cambia el valor de la referencia
    thetaRef=thetaRef*(-1); 
    kk=0;
 end
 TL=tL(ii);
 ref(ii)=thetaRef;
 estado=[ia(ii); omega(ii); theta(ii)];
 %Ley de control con ganancia de prealimentacion
 u(ii)=-K*estado+Gj*ref(ii);
 u(ii)=Max_U*tanh(u(ii)/Max_U) ; %Satura la tension de control.
 %Integracion Euler
 w_pp =(-w_p(ii)*(Ra*J+Laa*B)-omega(ii)*(Ra*B+Ki*Km)+u(ii)*Ki)/(J*Laa);
 ia_p=(1/Laa)*(-Ra*ia(ii)-Km*omega(ii)+u(ii));
 w_p(ii+1)=w_p(ii)+h*w_pp-(1/J)*TL;
 ia(ii+1)=ia(ii)+h*ia_p;
 omega(ii+1)=omega(ii)+h*w_p(ii);
 theta(ii+1)=theta(ii)+h*omega(ii);
 y_sal(ii)=Mat_C*estado;
 ii=ii+1;
end

%% Gráficos
fz=15;
t=0:h:tF;

figure(1);
legends_c = ["Obtenida", "Referencia"];
subplot(3,1,1);
plot(t,theta,'b','LineWidth', 1.2);
hold on;
plot(t,ref,'--r','LineWidth', 1.2);
grid on
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
title('Posicion Angular $\theta_t$ [rad]', 'Interpreter','latex','FontSize', fz);

subplot(3,1,2);
plot(t,ia,'r','LineWidth', 1.2);
title('Corriente $i_a$ [A]', 'Interpreter','latex','FontSize', fz);
grid on;

subplot(3,1,3);
plot(t,u,'m','LineWidth', 1.2);
title('Accion de control $u_t$ [V]', 'Interpreter','latex','FontSize', fz);
xlabel('Tiempo [seg]', 'Interpreter','latex','FontSize', fz-2); 
grid on;
ylim([-25 25])

set(gcf,'Color', 'w');
%%

legends_c = ["Obtenida", "Referencia"];
figure(2)
subplot(2,1,1);
title('Comparacion de Posicion Angular', 'Interpreter','latex','FontSize', fz);
plot(t,theta,'b','LineWidth', 1.2);
grid on;hold on; 
plot(t,ref,'--r','LineWidth', 1.4);
title('$\theta_t = \pi/2$', 'Interpreter','latex','FontSize', fz-1);
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
ylim([1.569 1.572])
xlim([0 0.6])
subplot(2,1,2);
plot(t,theta,'b','LineWidth', 1.2);
grid on;hold on; 
plot(t,ref,'--r','LineWidth', 1.4);
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
title('$\theta_t = -\pi/2$', 'Interpreter','latex','FontSize', fz-1);
xlabel('Tiempo [seg]', 'Interpreter','latex','FontSize', fz-2); 
ylim([-1.573 -1.569])
xlim([5 5.6])
set(gcf,'Color', 'w');
