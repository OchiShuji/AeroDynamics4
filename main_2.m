%/*************************************************************************
%* File Name     : main_2.m
%* Code Title    : 課題2:ゴドノフ法による線形双曲型連立偏微分方程式の解法
%* Programmer    :
%* Affiliation   : 
%* Creation Date : 2021/03/26
%* Language      : Matlab
%* Version       : 1.0.0
%**************************************************************************
%* NOTE:
%*  
%*************************************************************************/

clear all; close all;  %全ての変数とFigureを削除

%/*************************************************************************
tic()
A = [3.0,1.0;2.0,4.0];
[R,L] = eig(A);
l2 = L(2,2);
cfl = 0.5;
N = 50;
dx = 0.02;
dt = cfl * dx / l2;
travel = (dt*N).*L;

mx = 52;
x = zeros(1,mx);
xmin = -0.2;

U = zeros(2,mx);
W = zeros(2,mx);
U_init = zeros(2,mx);
W_init = zeros(2,mx);
U_exac = zeros(2,mx);
W_exac = zeros(2,mx);


%set grid
for i = 1:mx
    x(i) = xmin + dx*(i-2);
end

i_x0 = -xmin/dx+2;

%set initial condition
for i = 1:mx %境界外は0次補間
    if x(i) <= 0
        U_init(1,i) = 1;
        U_init(2,i) = 1;
    elseif x(i) > 0
        U_init(1,i) = 2;
        U_init(2,i) = -1;
    end
    W_init(:,i) = R\U_init(:,i);
end


%exact solution
for i = 2:mx-1
    if x(i) <= travel(1,1)
        W_exac(:,i) = W_init(:,i_x0-1);
    end
    if x(i) > travel(1,1) && x(i) <= travel(2,2)
        W_exac(1,i) = W_init(1,i_x0+1);
        W_exac(2,i) = W_init(2,i_x0-1);
    end
    if x(i) > travel(2,2) 
        W_exac(:,i) = W_init(:,i_x0+1);
    end
    U_exac(:,i) = R*W_exac(:,i);
end

F = zeros(2,mx);
F_bound = zeros(2,mx-1); 

for n = 1:N
    %F=AU
    for i = 1:mx
        if n == 1
            F(:,i) = A*U_init(:,i);
        elseif n > 1
            F(:,i) = A*U(:,i);
        end
    end
    %numerical flux(1st-order up-wind scheem)
    for i = 1:mx-1
        if n == 1
            F_bound(:,i) = 0.5.*(F(:,i)+F(:,i+1)-R*abs(L)*inv(R)*(U_init(:,i+1)-U_init(:,i)));
        elseif n > 1
            F_bound(:,i) = 0.5.*(F(:,i)+F(:,i+1)-R*abs(L)*inv(R)*(U(:,i+1)-U(:,i)));
        end 
    end
    %time integral(Godunov)
    for i = 2:mx-1
        if n == 1
            U(:,i) = U_init(:,i) - (dt/dx).*(F_bound(:,i)-F_bound(:,i-1));
        elseif n > 1
            U(:,i) = U(:,i) - (dt/dx).*(F_bound(:,i)-F_bound(:,i-1));
        end
    end
    %set boundary condition
    U(:,1) = U(:,2);
    U(:,mx) = U(:,mx-1);
end
toc()

figure(1);
plot(x(2:mx-1),U_init(1,2:mx-1),x(2:mx-1),U_init(2,2:mx-1),x(2:mx-1),U(1,2:mx-1),x(2:mx-1),U(2,2:mx-1),x(2:mx-1),U_exac(1,2:mx-1),x(2:mx-1),U_exac(2,2:mx-1));
xlim([-0.2 0.8]);
legend("u(initial)","v(initial)","u","v","u(exact)","v(exact)");
title("CFL="+cfl+", after "+N+" steps(t="+dt*N+")");


       
        


    