% Control con observador
% Control con Observador
clear all; clc; close all
%% Parámetros del sistema
Laa=4.994e-4; J=2.408e-9;Ra=20;B=0;Ki=9.647e-3;Km=60.529e-3;

%Paso, tiempo de simulación, pasos:
h=5e-7 ; tF=10; pasos=(tF/h);

%Definición e inicialización de variables
ia_p=0; w_p=0; w_pp=0; theta_p=0; 
ia_p_o=0; w_p_o=0; w_pp_o=0; theta_p_o=0; 
t=0:h:tF; ia=0:h:tF; omega=0:h:tF; w_p=0:h:tF; theta=0:h:tF;
ia_o=0:h:tF; omega_o=0:h:tF; w_p_o=0:h:tF; theta_o=0:h:tF;
u=linspace(0,0,pasos+1); ref=linspace(0,0,pasos+1); u_o=linspace(0,0,pasos+1);
thetaRef=(pi/2); %Referencia de ángulo del motor
Max_U=24; %Tension nominal del motor

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

tc=5; %Tiempo de cambio de la referencia (5s).
est=0; %Estado de la referencia.

ii=1; kk=0; 
%Condiciones iniciales
ia(1)=0; omega(1)=0; theta(1)=0; wp(1)=0; u(1)=0;
ia_o(1)=0; omega_o(1)=0; theta_o(1)=0; wp_o(1)=0; u_o(1)=0;
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
% mui(1)=-4e5; mui(2)=-5e4; mui(3)=-6e4; % Sin observador;
mui(1)=-7e4; mui(2)=-9e4; mui(3)=-25e3; 
alfa_i=poly(mui);
%Controlador
K=(fliplr(alfa_i(2:end)-c_ai(2:end)))*(inv(Mat_T));
disp('Polos de lazo cerrado controlador');
eig(Mat_A-(Mat_B*K)) %Polos de lazo cerrado
Gj=-inv(Mat_C*inv(Mat_A-(Mat_B*K))*Mat_B); %Ganancia de prealimentacion

%% Observador
Mat_A_O=Mat_A';
Mat_B_O=Mat_C';
Mat_C_O=Mat_B';

%Matríz Controlabilidad del observador:
Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O]; 
alfa_i_O=poly(26*mui); 
Mat_T_O=Mat_M_Dual*Mat_W;
Ko=(fliplr(alfa_i_O(2:end)-c_ai(2:end))*inv(Mat_T_O))';
disp('Polos de lazo cerrado observador');
eig(Mat_A-Ko*Mat_C) 

%Condicion inicial observador
x_hat=[0 ;0 ;0]; 

while(ii<(pasos+1))
    kk=kk+h;
    if(kk>tc) % Cada 5 s cambia el valor de la referencia
        thetaRef=thetaRef*(-1); 
        kk=0;
    end
    TL=tL(ii)/10;
    ref(ii)=thetaRef;
    estado_o=[ia_o(ii); omega_o(ii); theta_o(ii)];
    estado=[ia(ii); omega(ii); theta(ii)];
    
    %Ley de control
    u_o(ii)=-K*x_hat+Gj*ref(ii);        %Con observador
    u_o(ii)=Max_U*tanh(u_o(ii)/Max_U) ; %Satura la acción de control.
    
    u(ii)=-K*estado+Gj*ref(ii);         %Sin observador
    u(ii)=Max_U*tanh(u(ii)/Max_U) ;     %Satura la tension de control.
    
    %Integracion Euler sin observador:
    w_pp =(-w_p(ii)*(Ra*J+Laa*B)-omega(ii)*(Ra*B+Ki*Km)+u(ii)*Ki)/(J*Laa);
    ia_p=(1/Laa)*(-Ra*ia(ii)-Km*omega(ii)+u(ii));
    w_p(ii+1)=w_p(ii)+h*w_pp-(1/J)*TL;
    ia(ii+1)=ia(ii)+h*ia_p;
    omega(ii+1)=omega(ii)+h*w_p(ii);
    theta(ii+1)=theta(ii)+h*omega(ii);
 
    %Integracion Euler con observador:
    w_pp_o =(-w_p_o(ii)*(Ra*J+Laa*B)-omega_o(ii)*(Ra*B+Ki*Km)+u_o(ii)*Ki)/(J*Laa);
    ia_p_o=(1/Laa)*(-Ra*ia_o(ii)-Km*omega_o(ii)+u_o(ii));
    w_p_o(ii+1)=w_p_o(ii)+h*w_pp_o-(1/J)*TL;
    ia_o(ii+1)=ia_o(ii)+h*ia_p_o;
    omega_o(ii+1)=omega_o(ii)+h*w_p_o(ii);
    theta_o(ii+1)=theta_o(ii)+h*omega_o(ii);
    
    % Observador
    y_sal(ii)=Mat_C*estado_o;
    y_sal_O(ii)=Mat_C*x_hat;
    x_hatp=Mat_A*x_hat+Mat_B*u_o(ii)+Ko*(y_sal(ii)-y_sal_O(ii));
    x_hat=x_hat+h*x_hatp;
 
    ii=ii+1;
end

%% Gráficos
fz=15;
t=0:h:tF;

figure(1);
legends_c = ["Obtenida", "Referencia"];
subplot(3,1,1);
plot(t,theta_o,'b','LineWidth', 1.2);
hold on;
plot(t,ref,'--r','LineWidth', 1.2);
grid on
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
title('Posicion Angular $\theta_t$ [rad]', 'Interpreter','latex','FontSize', fz);

subplot(3,1,2);
plot(t,ia_o,'r','LineWidth', 1.2);
title('Corriente $i_a$ [A]', 'Interpreter','latex','FontSize', fz);
xlabel('Tiempo en Seg.', 'Interpreter','latex','FontSize', fz-2); 
grid on;

subplot(3,1,3);
plot(t,u_o,'m','LineWidth', 1.2);
title('Accion de control $u_t$ [V]', 'Interpreter','latex','FontSize', fz);
xlabel('Tiempo en Seg.', 'Interpreter','latex','FontSize', fz-2); 
grid on;
ylim([-25 25])

set(gcf,'Color', 'w');
%%

legends_c = ["Obtenida sin Observador","Obtenida con Observador","Referencia"];

figure(2)
title('Comparacion de Posicion Angular', 'Interpreter','latex','FontSize', fz);
subplot(2,1,1);
plot(t,theta,'r','LineWidth', 1.2);
hold on; 
plot(t,theta_o,'b','LineWidth', 1.2);
plot(t,ref,'--k','LineWidth', 1);
grid on;
title('$\theta_t = \pi/2$', 'Interpreter','latex','FontSize', fz-1);
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
ylim([1.56 1.575])
xlim([0 0.6])
subplot(2,1,2);
plot(t,theta,'r','LineWidth', 1.2);
hold on;
plot(t,theta_o,'b','LineWidth', 1.2);
plot(t,ref,'--k','LineWidth', 1);
grid on; 
title('$\theta_t = -\pi/2$', 'Interpreter','latex','FontSize', fz-1);
legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-2);
xlabel('Tiempo en Seg.', 'Interpreter','latex','FontSize', fz-2); 
ylim([-1.58 -1.565])
xlim([5 5.6])
set(gcf,'Color', 'w');
