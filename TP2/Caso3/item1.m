% Sin observador
clc; clear all;
%Parámetros del sistema
m=0.1; Fricc=0.1; long=1.6; g=9.8; M=1.5;
%Tiempo de muestreo, paso de integracion, tiempo de simulación:
Ts=0.001 ; VecesEuler=10; h=Ts/VecesEuler; T=50 ; KMAX=T/Ts; 
%Definición, Inicialización de variables
p_pp=0; fi_pp=0; t=0:h:T;
ii=1; flag=0;
%Versión linealizada en el equilibrio inestable:
Mat_Ac=[0 1 0 0; 0 -Fricc/M -m*g/M 0; 0 0 0 1;0 -Fricc/(long*M) -g*(m+M)/(long*M) 0]
Mat_Bc=[0; 1/M; 0; 1/(long*M)]
Mat_C =[1 0 0 0; %La salida son la posicion y el angulo fi
 0 0 1 0];
Cref=Mat_C(1,:) ; %Extrae la primer fila de C.
Mat_D = [0 ; 0];
auto_valc=eig(Mat_Ac)
%% DISCRETIZACIÓN
sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,Mat_D);
sys_d=c2d(sys_c,Ts,'zoh'); 
Mat_A=sys_d.a; 
Mat_B=sys_d.b;
%Polos de lazo abierto
auto_val=eig(Mat_A)
%% SISTEMA AMPLIADO
Mat_Aa = [Mat_A ,zeros(4,1) ; -Cref*Mat_A ,1]
Mat_Ba = [Mat_B ; -Cref*Mat_B]
Mat_Ma=[Mat_Ba Mat_Aa*Mat_Ba Mat_Aa^2*Mat_Ba Mat_Aa^3*Mat_Ba Mat_Aa^4*Mat_Ba];%Matriz Controlabilidad
rango=rank(Mat_Ma)
%% CONTROLADOR DLQR
d=[1 0.001 1 1 .0000001];
d=[12000 8000 1e-2 5000 .0008];
d=[80 1000 1 1 .000005];
Q=diag(d); 
R=1000; 
[Ka,P,E] = dlqr(Mat_Aa,Mat_Ba,Q,R); %E:Vector de valores propios
 %P:Matriz simetrica definida positiva
disp('Polos lazo cerrado del controlador ampliado: '); 
eig(Mat_Aa-Mat_Ba*Ka) %Polos de lazo cerrado
K=Ka(1:end-1) %Controlador
KI=-Ka(end) %Constante de ganancia integral
disp('Polos lazo cerrado del controlador: ');
aut_controlador=abs(eig(Mat_Aa-Mat_Ba*Ka))
%% Condiciones Iniciales
ref_pos=10; % Referencia delta
p(1)=0; p_p(1)=0; fi_p(1)=0; u(1)=0; ve(1)=0;
fi(1)=pi; color='r';
xo=[0;0;pi;0]; %Condicion incial equilibrio estable
%Vector de estados
x=[p(1); p_p(1); fi(1); fi_p(1)];
%% SIMULACION
for ki=1:KMAX
    %Ley de control
    y_sal=Mat_C*x;
    ve(ki+1)=ve(ki)+ref_pos-y_sal(1);
    u1(ki)=-K*x+KI*ve(ki+1); color='r';
    % Integracion de Euler
    for kii=1:Ts/h
        u(ii)=u1(ki);
        p_pp=(1/(M+m))*(u(ii)-m*long*fi_pp*cos(fi(ii))+m*long*(fi_p(ii)^2)*sin(fi(ii))-Fricc*p_p(ii));
        fi_pp=(1/long)*(g*sin(fi(ii))-p_pp*cos(fi(ii)));
        p_p(ii+1)=p_p(ii)+h*p_pp;
        p(ii+1)=p(ii)+h*p_p(ii);
        fi_p(ii+1)=fi_p(ii)+h*fi_pp;
        fi(ii+1)=fi(ii)+h*fi_p(ii);
        %Cambio de referencia de desplazamiento y de masa
        epsilon=(abs(p(ii)-10))^2+(abs(p_p(ii)-0))^2;
        if(epsilon<0.001)
            if(flag==0)
                ref_pos=0; %Cambia la referencia de desplazamiento a 0 metros
                m=m*10; %Aumenta la masa 10 veces
                flag=1; 
            end
        end
        ii=ii+1; 
    end
    x=[p(ii-1); p_p(ii-1); fi(ii-1) ; fi_p(ii-1) ];
end
u(ii)=u1(ki);

%% Gráficos
fz=15;

figure(1);

subplot(3,2,1);
plot(t,fi_p,color,'Linewidth',1.2);
grid on; 
title('Velocidad angular [rad/s]', 'Interpreter','latex','FontSize', fz);
hold on;

subplot(3,2,2);plot(t,fi,color,'Linewidth',1.2);
grid on;
title('$\phi$ , Angulo [rad]', 'Interpreter','latex','FontSize', fz);
hold on;
plot(t,pi*ones(size(t)),'LineStyle','--','color','k');
hold on;

subplot(3,2,3); plot(t,p,color,'Linewidth',1.2);
grid on;
title('$\delta$ , Posicion carro [m]', 'Interpreter','latex','FontSize', fz);
ylim([-1 12]);hold on;

subplot(3,2,4);
plot(t,p_p,color,'Linewidth',1.2);
grid on;
title('Velocidad carro [m/s]', 'Interpreter','latex','FontSize', fz);
hold on;

subplot(3,1,3);
plot(t,u,color,'Linewidth',1.2);
grid on;title('u , Accion de control', 'Interpreter','latex','FontSize', fz);
xlabel('Tiempo [s]', 'Interpreter','latex','FontSize', fz-2);
hold on;

set(gcf,'Color', 'w');

figure(2);hold on;

subplot(2,1,1);
plot(fi,fi_p,color,'Linewidth',1.2);
grid on;
xlabel('Angulo', 'Interpreter','latex','FontSize', fz-2);
ylabel('Velocidad angular', 'Interpreter','latex','FontSize', fz-2);
hold on;

subplot(2,1,2);
plot(p,p_p,color,'Linewidth',1.2);
grid on;
xlabel('Posicion carro', 'Interpreter','latex','FontSize', fz-2);
ylabel('Velocidad carro', 'Interpreter','latex','FontSize', fz-2);
hold on;

set(gcf,'Color', 'w');
