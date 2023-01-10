% this Program for : 1-Dim., unsteady State,un-lumped system, without inner heat generation.
clear;clc;clf;
% ==== Input Data ====

M = 10 ;              % no. of points

L = .12;          %  wall thickness [m]
k = 10;              % thermal conductivity [w/m.k]
alpha = 2.5*10^(-6);   % thermal conductivity [m^2/s]
q = 10^4;            % constant heat flux[w/m^2]
h = 700;              % [w/m^2.k]
T_inf = 300;         % [c]
Ti = 20 ;             %initial temperature at time=0 [c]
t = 30*60;           % total time [sec.]

% ==== manual calculation ====

dx = L/(M-1);
tao_min = 1/(2+(2*dx*h/k));
dt = ceil(((dx^2)/alpha)*tao_min); % [sec.] %ceil funtion to rounding decimal fractions
tao = (alpha/(dx^2))*dt;
time_max=1800;
I = ceil(time_max/dt);
x = (0:dx:L).*100 ;     % [cm]
% ==== initial state ====

T = zeros(I,M);
T(1,:)= Ti;



% ==== Numerical ====


for i=1:I;
for m=1:M;
    if m==1
       T(i+1,m) = (1-2*tao-((2*dx*h*tao)/k))*T(i,m)+((2*dx*h*tao)/k)*T_inf+2*tao*T(i,m+1) ; 
    elseif m==M
        T(i+1,m) = 2*tao*T(i,m-1)+(1-2*tao)*T(i,m)-(q*dx*2*tao)/k;
    else
        T(i+1,m) = (1-2*tao)*T(i,m)+tao*T(i,m-1)+tao*T(i,m+1);
    end
end

end



                       % ==== Outputs ====

% ==== A ====
y5 = T(ceil((5*60)/dt),1:M);
y10 = T(ceil((10*60)/dt),1:M);
y15 = T(ceil((15*60)/dt),1:M);
y20 = T(ceil((20*60)/dt),1:M);
y25 = T(ceil((25*60)/dt),1:M);
y30 = T(ceil(I),1:M);

Min_temp = T(:,M);
time1=(0:(t/I):t)/60;

figure(1)
plot(x,y5,'--b','linewidth',1)
hold on
plot(x,y5,'--k','linewidth',1)
plot(x,y10,'--b','linewidth',1)
plot(x,y15,'--m','linewidth',1)
plot(x,y20,'--g','linewidth',1)
plot(x,y25,'--b','linewidth',1)
plot(x,y30,'-r','linewidth',1)
hold off


grid on;grid minor;box on
xlabel('x  (cm)','fontsize',20)
ylabel('Temp.  (c)','fontsize',20)
title('temperature distribution','fontsize',20)
legend('temp at 5min' , 'temp at 10min','temp at 15min','temp at 20min','temp at 25min','temp at 30min')
saveas(gcf, 'no2');


% ==== B ====
Min_temp = T(:,M);
time1=(0:(t/I):t)/60;

min = y30==min(y30);
fprintf('after %0.0f min. \n',t/60)
fprintf('min. temp. = %0.4f [c]\n',y30(min))
fprintf('at X = %0.4f [cm]\n',x(min))


figure(2)
plot(time1,Min_temp,'-b','linewidth',1)
grid on;grid minor;box on
xlabel('time  (sec)','fontsize',20)
ylabel('min_Temp.  (c)','fontsize',20)
title('min Temp. vs time','fontsize',20)
saveas(gcf, 'no2');
 

% ==== C ====
for i=I ,
for  m=1:M   
    if m==1
        Q_1 = h*(T_inf-T(i,m))- (k*(T(i,m)-T(i,m+1)))/dx ;       % where Q@1  is  stored energy at external node 1  
    elseif m==M     
        Q_M = (k*(T(i,m-1)-T(i,m)))/dx - q ;                     % where Q@M  is  stored energy at external node M
    end
end
end
        
Q_total = Q_1 + Q_M;                  % where Q_total is  total stored energy after 30 min [J] 

fprintf('the total amount of heat stored in the plate after 30 min =   %0.4f [KJ]\n',Q_total/10^3)
 
