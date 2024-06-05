%Zona muerta
clc; clear all;
%Parámetros del sistema
m=0.1; Fricc=0.1; long=1.6; g=9.8; M=1.5; m_o=m;
%Tiempo de muestreo, paso de integracion, tiempo de simulación:
Ts=0.001 ; VecesEuler=10; h=Ts/VecesEuler; T=50 ; KMAX=T/Ts; 
%Definición, Inicialización de variables
t=0:h:T; p_pp_o=0; fi_pp_o=0;
jj=1; flag_o=0; 
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
%% OBSERVADOR
Mat_A_O=Mat_A';
Mat_B_O=Mat_C';
Mat_C_O=Mat_B';

% CÁLCULO Ko POR DLQR
d=[1 100 0.01 1000];
Qo=diag(d); 
Ro=diag([1 1e6]); 
[Ko,Po,Eo] = dlqr(Mat_A_O,Mat_B_O,Qo,Ro); %E:Vector de valores propios
 %P:Matriz simetrica definida positiva
Kobs=Ko'
disp('Polos de lazo cerrado del observador: ');
p_observador=abs(eig(Mat_A-Kobs*Mat_C)) %Verifico que todos los polos estén en el semiplano izquierdo
%% Condiciones Iniciales

ref_pos_o=10; % Referencia delta
p_o(1)=0; p_p_o(1)=0; fi_p_o(1)=0; u_o(1)=0; ve_o(1)=0;
fi_o(1)=pi;
xo=[0;0;pi;0]; %Condicion incial equilibrio estable
%Vector de estados
x_o=[p_o(1); p_p_o(1); fi_o(1); fi_p_o(1)];
x_hat=[0;0;0;0]; %Inicializo el Observador
%% SIMULACION
for ki=1:KMAX
    %Ley de control
    y_sal_o=Mat_C*x_o;
    ve_o(ki+1)=ve_o(ki)+ref_pos_o-y_sal_o(1);
    u1_o(ki)=-K*x_hat+KI*ve_o(ki+1); %Con observador de estados
    
    % No linealidad
    ZM=0.01;
    if(abs(u1_o(ki))<ZM)
        u1_o(ki)=0;
    else
        u1_o(ki)=sign(u1_o(ki))*(abs(u1_o(ki))-ZM);
    end
    
    % Integracion de Euler para sistema con observador:
    for kii=1:Ts/h
        u_o(jj)=u1_o(ki);
        p_pp_o=(1/(M+m_o))*(u_o(jj)-m_o*long*fi_pp_o*cos(fi_o(jj))+m_o*long*(fi_p_o(jj)^2)*sin(fi_o(jj))-Fricc*p_p_o(jj));
        fi_pp_o=(1/long)*(g*sin(fi_o(jj))-p_pp_o*cos(fi_o(jj)));
        p_p_o(jj+1)=p_p_o(jj)+h*p_pp_o;
        p_o(jj+1)=p_o(jj)+h*p_p_o(jj);
        fi_p_o(jj+1)=fi_p_o(jj)+h*fi_pp_o;
        fi_o(jj+1)=fi_o(jj)+h*fi_p_o(jj);
        %Cambio de referencia de desplazamiento y de masa
        epsilon=(abs(p_o(jj)-10))^2+(abs(p_p_o(jj)-0))^2;
        if(epsilon<0.001)
            if(flag_o==0)
                ref_pos_o=0; %Cambia la referencia de desplazamiento a 0 metros
                m_o=m*10; %Aumenta la masa 10 veces
                flag_o=1; 
            end
        end
        jj=jj+1; 
    end
    
    x_o=[p_o(jj-1); p_p_o(jj-1); fi_o(jj-1) ; fi_p_o(jj-1) ];
    
    y_hat=Mat_C*x_hat;
    x_hat=Mat_A*x_hat+Mat_B*u1_o(ki)+Kobs*(y_sal_o-y_hat); %Aca se usa y. 
end
u_o(jj)=u1_o(ki);

%% Gráficos
fz=15;


color = [1 0 0];
color = [0 0 1];
color = [1 0 1];
color = [0.4660 0.6740 0.1880];
color = [0.4940 0.1840 0.5560];

figure(1);
legends_c = ["ZM=0.5","ZM=0.3","ZM=0.2","ZM=0.05","ZM=0.01"];
legends_c2 = ["ZM=0.5","ZM=0.3","ZM=0.2","ZM=0.05","ZM=0.01","Referencia"];

subplot(3,2,1);
hold on;
plot(t,fi_p_o,'Color',color,'Linewidth',1.2);
grid on; 
title('Velocidad angular [rad/s]', 'Interpreter','latex','FontSize', fz);
%legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-5);

subplot(3,2,2);
hold on;
plot(t,fi_o,'Color',color,'Linewidth',1.2);
grid on;
title('$\phi$ , Angulo [rad]', 'Interpreter','latex','FontSize', fz);
%plot(t,pi*ones(size(t)),'LineStyle','--','color','k');
%legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-5);
hold on;

subplot(3,2,3);
hold on;
plot(t,p_o,'Color',color,'Linewidth',1.2);
grid on;
title('$\delta$ , Posicion carro [m]', 'Interpreter','latex','FontSize', fz);
%legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-5);
ylim([-1 12]);

subplot(3,2,4);
hold on;
plot(t,p_p_o,'Color',color,'Linewidth',1.2);
grid on;
title('Velocidad carro [m/s]', 'Interpreter','latex','FontSize', fz);
%legend(legends_c,'Location','northeast','Interpreter','latex','FontSize', fz-5);

subplot(3,1,3);
hold on;
plot(t,u_o,'Color',color,'Linewidth',1.2);
grid on;title('u , Accion de control', 'Interpreter','latex','FontSize', fz);
xlabel('Tiempo [s]', 'Interpreter','latex','FontSize', fz-2);
hold on;
%legend(legends_c,'Location','southeast','Interpreter','latex','FontSize', fz-5);

set(gcf,'Color', 'w');
