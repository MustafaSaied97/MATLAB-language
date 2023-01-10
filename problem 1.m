% this Program for : 1-Dim., Steady State, with inner heat generation. 
clc;clf;clear;

% ==== Input Data ====
M = 100;        % no. of points

L = 0.2;     % [m]
dx = L/(M-1); % step
x = 0:dx:L;
k = 15;      %thermal conductivity [w/m.k]
g_o = 10^5;     
a = 3;
g = g_o*exp(-a*x/L);    %heat generation function [w/m^3]
h = 30;       % heat transfer coefficient  [w/m^2.k]
T_inf = 20+273;  % [k]
q_rad = 200;      % [w/m^2]

% ==== Numerical ====

T = zeros(1,M)+100;     % initail guess

error = 1;
num = 0;

while error > 0.0001;
    num = num+1;
    for m = 1:M
        if m==1
            temp = T(m);
            T(m) = T(m+1)+(dx^2).*g(m)./(2*k);
            err(m) = abs(temp-T(m));
        else if m==M
                temp = T(m);
                T(m) = ((k/dx)*T(m-1)+h*T_inf-q_rad+dx*g(m)/2)/(k/dx+h);
                err(m) = abs(temp-T(m));
        else 
            temp = T(m);
            T(m) = 0.5*(T(m-1)+T(m+1)+(dx^2)*g(m)/k);
            err(m) = abs(temp-T(m));
            end;end;end
    error = max(err);
end

% ==== Outputs ====


plot(x,T,'--b','linewidth',1);


grid on;grid minor;box on
xlabel('x  (m)','fontsize',20)
ylabel('Temp.  (k)','fontsize',20)
title('No. 1 temp distribution','fontsize',20)


fprintf('Error = %0.4e\n',error)
fprintf('no. of Iteration = %0.0f\n',num)


t_max = T==max(T);
fprintf('max. temp. = %0.4f [k]\n',T(t_max))
fprintf('at X = %0.4f [m]\n',x(t_max))

