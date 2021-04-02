%/*************************************************************************
%* File Name     : main_3.m
%* Code Title    : 課題3:衝撃波管問題の数値解法
%* Programmer    :
%* Affiliation   : 
%* Creation Date : 2021/03/30
%* Language      : Matlab
%* Version       : 1.0.0
%**************************************************************************
%* NOTE:
%*  
%*************************************************************************/

clear all; close all;  %全ての変数とFigureを削除

%/*************************************************************************

%set grid
mx = 400;
x  = zeros(1,mx+5);
xmin = -2.0;
xmax = 2.0;
dx = (xmax-xmin)/mx;
for i = 1:mx+5
    x(i) = xmin + dx*(i-3);
end

%Set parameters
p_atm = 1.013e+5; %1[atm] = 1.013e+5[Pa]
eps = 1.0e-6;
gamma = 1.4; %Cp/Cv
R = 8.314e+3/28.96; %gas constant[m^2/s^2/K]
C_v = R/(gamma-1);

%initial value input
u_l = 0.0; %[m/s]
u_r = 0.0; %[m/s]
t_l = 400; %[K]
t_r = 300; %[K]
p_l = 5*p_atm; %[Pa]
p_r = 1*p_atm; %[Pa]
rho_r = p_r / (R*t_r); %[kg/m^3]
rho_l = p_l / (R*t_l); %[kg/m^3]

cfl = 0.25;

ecp = 0.125; %entropy correction parameter

u_ref = sqrt(gamma*R*t_l); %reference speed

dt = cfl*dx/sqrt(gamma*R*max(t_l,t_r));
disp("dt="+dt*1000+"[ms]");
nlast = 400; %final time step
time = nlast*dt;
disp("time="+time*1000+"[ms]");

u = zeros(1,mx+5,nlast);
rho = zeros(1,mx+5,nlast);
e = zeros(1,mx+5,nlast);
t = zeros(1,mx+5,nlast);
p = zeros(1,mx+5,nlast);
Q = zeros(3,mx+5); %conservative variables
E = zeros(3,mx+5); %flux

%set initial condition
%領域外もset,0次補間
for i = 1:mx+5
    if x(i) < 0
        rho(1,i,1) = rho_l;
        u(1,i,1) = u_l;
        p(1,i,1) = p_l;
        t(1,i,1) = t_l;
    elseif x(i) >= 0
        rho(1,i,1) = rho_r;
        u(1,i,1) = u_r;
        p(1,i,1) = p_r;
        t(1,i,1) = t_r;
    end
    e(1,i,1) = p(1,i,1)/(rho(1,i,1)*(gamma-1));
    Q(1,i) = rho(1,i,1);
    Q(2,i) = rho(1,i,1)*u(1,i,1);
    Q(3,i) = rho(1,i,1)*(e(1,i,1)+0.5*u(1,i,1)^2); %Et:Total Energy
end

%time marching
Q_old = zeros(3,mx+5);
for n = 1:nlast-1
    Q_old = Q;
    [E] = calc_flux(Q,mx,gamma,E,ecp,dt,dx);
    for i = 3:mx+3
        Q(:,i) = Q_old(:,i) - (dt/dx).*(E(:,i)-E(:,i-1));
        rho(1,i,n+1) = Q(1,i);
        u(1,i,n+1) = Q(2,i)/rho(1,i,n+1);
        p(1,i,n+1) = (gamma-1.0)*(Q(3,i)-0.5*rho(1,i,n+1)*u(1,i,n+1)^2);
        e(1,i,n+1) = p(1,i,n+1)/((gamma-1.0)*rho(1,i,n+1));
        t(1,i,n+1) = e(1,i,n+1)/C_v;
    end
    disp(n+" steps calculated, t="+dt*n);
end

%% 結果表示
figure(1);
plot(x(3:mx+3),p(1,3:mx+3,nlast));
legend("p");
title("CFL="+cfl+", t="+time);
xlabel("x[m]");
ylabel("pressure[Pa]");

figure(2);
plot(x(3:mx+3),rho(1,3:mx+3,nlast));
legend("rho");
title("CFL="+cfl+", t="+time);
xlabel("x[m]");
ylabel("density[kg/m^3]");

figure(3);
plot(x(3:mx+3),u(1,3:mx+3,nlast));
legend("u");
title("CFL="+cfl+", t="+time);
xlabel("x[m]");
ylabel("velocity[m/s]");

t_step = zeros(1,nlast);
for n = 1:nlast
    t_step(n) = 1000*dt*(n-1);
end

figure(4);
[T,X] = meshgrid(t_step,x(3:mx+3));
P = reshape(p,[mx+5,nlast]);
contour(X,T,P(3:mx+3,1:nlast),20);
title("Pressure Wave Diagram(CFL="+cfl);
xlabel("x[m]");
ylabel("t[ms]");

figure(5);
[T,X] = meshgrid(t_step,x(3:mx+3));
U = reshape(u,[mx+5,nlast]);
contour(X,T,U(3:mx+3,1:nlast),20);
title("Velocity Wave Diagram(CFL="+cfl);
xlabel("x[m]");
ylabel("t[ms]");

figure(6);
[T,X] = meshgrid(t_step,x(3:mx+3));
TEMP = reshape(t,[mx+5,nlast]);
contour(X,T,TEMP(3:mx+3,1:nlast),20);
title("Temperature Wave Diagram(CFL="+cfl);
xlabel("x[m]");
ylabel("t[ms]");

%% fluxの計算
function [E] = calc_flux(Q,mx,gamma,E,ecp,dt,dx)
lambda = zeros(3,mx+5);
R = zeros(3,3,mx+5);
R_inv =zeros(3,3,mx+5);
alpha = zeros(3,mx+5);

for i = 1:mx+4
    % i:left,i+1:right
    dq = Q(:,i+1)-Q(:,i);
    
    rho_lt = Q(1,i);
    u_lt = Q(2,i)/rho_lt;
    p_lt = (gamma-1.0)*(Q(3,i)-0.5*rho_lt*u_lt^2);
    h_lt = (Q(3,i)+p_lt)/rho_lt; %total enthalpy(h=(Et+P)/rho)
    
    rho_rt = Q(1,i+1);
    u_rt = Q(2,i+1)/rho_lt;
    p_rt = (gamma-1.0)*(Q(3,i+1)-0.5*rho_rt*u_rt^2);
    h_rt = (Q(3,i+1)+p_rt)/rho_rt;
    
    %Roe average
    u_bar = (sqrt(rho_lt)*u_lt+sqrt(rho_rt)*u_rt)/(sqrt(rho_lt)+sqrt(rho_rt));
    h_bar = (sqrt(rho_lt)*h_lt+sqrt(rho_rt)*h_rt)/(sqrt(rho_lt)+sqrt(rho_rt));
    a_bar = sqrt(max((gamma-1.0)*(h_bar-0.5*u_bar^2),min(gamma*p_lt/rho_lt,gamma*p_rt/rho_rt)));
    
    lambda(:,i) = [u_bar-a_bar;
                   u_bar;
                   u_bar+a_bar];
    R(:,:,i) = [1.0,1.0,1.0;
                u_bar-a_bar,u_bar,u_bar+a_bar;
                h_bar-u_bar*a_bar,0.5*u_bar^2,h_bar+u_bar*a_bar];
            
    b2 = (gamma-1.0)/a_bar^2;
    b1 = b2*u_bar^2/2.0;
    
    R_inv(:,:,i) = [0.5*(b1+u_bar/a_bar),0.5*(-b2*u_bar-1/a_bar),0.5*b2;
                    1-b1,b2*u_bar,-b2;
                    0.5*(b1-u_bar/a_bar),0.5*(-b2*u_bar+1/a_bar),0.5*b2];
    %TVD scheme
    alpha(:,i) = R_inv(:,:,i)*dq;
end

for i = 2:mx+3
    phi = zeros(3,1);
    tvd = zeros(3,1);
    ecpx = ecp*max(abs(lambda(:,i)));
    %for lambda1 = u-c
    ff1 = (dt/dx)*lambda(1,i)^2;
    ff2 = abs(lambda(1,i));
    if ff2 < ecpx
        ff2 = (ff2^2+ecpx^2)*0.5/ecpx;  %entropy correction
    end
    sig = sign(alpha(1,i));
    q_limit = sig*max(0.0,min([sig*2.0*alpha(1,i-1),sig*2.0*alpha(1,i),sig*2.0*alpha(1,i+1),sig*0.5*(alpha(1,i-1)+alpha(1,i+1))]));
    phi(1,1) = -ff1*q_limit-ff2*(alpha(1,i)-q_limit);
    %for lambda2 = u
    ff1 = (dt/dx)*lambda(2,i)^2;
    ff2 = abs(lambda(2,i));
    if ff2 < ecpx
        ff2 = (ff2^2+ecpx^2)*0.5/ecpx;  %entropy correction
    end
    sig = sign(alpha(2,i));
    q_limit = sig*max(0.0,min([sig*2.0*alpha(2,i-1),sig*2.0*alpha(2,i),sig*2.0*alpha(2,i+1),sig*0.5*(alpha(2,i-1)+alpha(2,i+1))]));
    phi(2,1) = -ff1*q_limit-ff2*(alpha(2,i)-q_limit);
    
    %for lambda2 = u+c
    ff1 = (dt/dx)*lambda(3,i)^2;
    ff2 = abs(lambda(3,i));
    if ff2 < ecpx
        ff2 = (ff2^2+ecpx^2)*0.5/ecpx;  %entropy correction
    end
    sig = sign(alpha(3,i));
    q_limit = sig*max(0.0,min([sig*2.0*alpha(3,i-1),sig*2.0*alpha(3,i),sig*2.0*alpha(3,i+1),sig*0.5*(alpha(3,i-1)+alpha(3,i+1))]));
    phi(3,1) = -ff1*q_limit-ff2*(alpha(3,i)-q_limit);
    
    tvd(:,1) = R(:,:,i)*phi(:,1);
    
    %fluxの計算
    rho_lt = Q(1,i);
    u_lt = Q(2,i)/rho_lt;
    p_lt = (gamma-1.0)*(Q(3,i)-0.5*rho_lt*u_lt^2);
    
    rho_rt = Q(1,i+1);
    u_rt = Q(2,i+1)/rho_lt;
    p_rt = (gamma-1.0)*(Q(3,i+1)-0.5*rho_rt*u_rt^2);
    
    E_l = [rho_lt*u_lt;rho_lt*u_lt^2+p_lt;u_lt*(Q(3,i)+p_lt)];
    E_r = [rho_rt*u_rt;rho_rt*u_rt^2+p_rt;u_rt*(Q(3,i)+p_rt)];
    
    E(:,i) = 0.5*(E_r+E_l+tvd(:,1));
end
end





