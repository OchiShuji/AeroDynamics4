%/*************************************************************************
%* File Name     : main_1.m
%* Code Title    : 課題1:スカラー1次元線形双曲型偏微分方程式
%* Programmer    :
%* Affiliation   : 
%* Creation Date : 2020/06/13
%* Language      : Matlab
%* Version       : 1.0.0
%**************************************************************************
%* NOTE:
%*  
%*************************************************************************/

clear all; close all;  %全ての変数とFigureを削除

%/*************************************************************************
tic()
x = zeros(1,102);
yi = zeros(1,102); %initial state
ye = zeros(1,102); %exact solution
flux = zeros(1,102); %numerical flux
work = zeros(1,102); 
ys1 = zeros(1,102); %solution by minmod
ys2 = zeros(1,102); %solution by superbee
ys3 = zeros(1,102); %solution by Sym-TVD(minmod)
ys4 = zeros(1,102); %solution by Sym-TVD(superbee)
ys5 = zeros(1,102); %solution by 1st order upwind scheme
ys6 = zeros(1,102); %solution by Lax-Wendroff scheme

c = 1.0;
mx = 501;
nlast = 300;
cfl = 0.5;
dx = 0.002;
dt = cfl*dx;
travel = dt*c*nlast;

%set grid
for i = 1:mx+3
    x(i) = dx*(i-3);
end

%set init condition
for i = 1:mx+3
    if x(i) < 0.5
        yi(i) = 1.0;
    else
        yi(i) = 0.0;
    end
end

%set exact solution
for i = 1:mx+3
    if x(i) < 0.5+travel
        ye(i) = 1.0;
    else
        ye(i) = 0.0;
    end
end

%% TVD scheme
%solve scheme by minimod
TV1 = 0.0;
for i = 1:mx+3
    ys1(i) = yi(i);
end

for n = 1:nlast
    for i = 2:mx+1  
        dm = ys1(i)-ys1(i-1);
        d0 = ys1(i+1)-ys1(i);
        dp = ys1(i+2)-ys1(i+1);
        B = minmod(1,dm/d0);
        flux(i) = c*(ys1(i)+0.5*(1-cfl)*B*d0);
    end
    for i = 3:mx+1
        work(i) = ys1(i)-(dt/dx)*(flux(i)-flux(i-1));
    end
    for i = 3:mx+1
        ys1(i) = work(i);
    end
end
for i = 4:mx+2
    TV1 = TV1+abs(ys1(i)-ys1(i-1));
end

%solve scheme by superbee
TV2 = 0.0;
for i = 1:mx+3
    ys2(i) = yi(i);
end

for n = 1:nlast
    for i = 2:mx+1  
        dm = ys2(i)-ys2(i-1);
        d0 = ys2(i+1)-ys2(i);
        dp = ys2(i+2)-ys2(i+1);
        B = max(0.0,max(min(2*dm/d0,1),min(2,dm/d0)));
        flux(i) = c*(ys2(i)+0.5*(1-cfl)*B*d0);
    end
    for i = 3:mx+1
        work(i) = ys2(i)-(dt/dx)*(flux(i)-flux(i-1));
    end
    for i = 3:mx+1
        ys2(i) = work(i);
    end
end
for i = 4:mx+2
    TV2 = TV2+abs(ys2(i)-ys2(i-1));
end
%% Yee's symmetric TVD scheme
%solve scheme by minimod
ecp = 0.1; %entropy correction parameter
TV3 = 0.0;

for i = 1:mx+3
    ys3(i) = yi(i);
end

for n = 1:nlast
    for i = 2:mx+1  
        dm = ys3(i)-ys3(i-1);
        d0 = ys3(i+1)-ys3(i);
        dp = ys3(i+2)-ys3(i+1);
        q = minmod(dm,minmod(d0,dp));
        ac = abs(c);
        if ac < ecp
            %entropy correction
            ac = (c^2+ecp^2)*0.5/ecp;
        end
        f = -(dt*c^2/dx*q+ac*(d0-q));
        flux(i) = 0.5*(c*ys3(i)+c*ys3(i+1)+f);
    end
    for i = 3:mx+1
        work(i) = ys3(i)-(dt/dx)*(flux(i)-flux(i-1));
    end
    for i = 3:mx+1
        ys3(i) = work(i);
    end
end
for i = 4:mx+2
    TV3 = TV3+abs(ys3(i)-ys3(i-1));
end

%solve scheme by superbee
TV4 = 0.0;
for i = 1:mx+3
    ys4(i) = yi(i);
end

for n = 1:nlast
    for i = 2:mx+1  
        dm = ys4(i)-ys4(i-1);
        d0 = ys4(i+1)-ys4(i);
        dp = ys4(i+2)-ys4(i+1);
        s = sign(dm);
        sb1 = s*max(0.0,max(min(2*abs(dm),s*d0),min(abs(dm),2*s*d0)));
        s = sign(d0);
        sb2 = s*max(0.0,max(min(2*abs(d0),s*dp),min(abs(d0),2*s*dp)));
        q = sb1+sb2-d0;
        ac = abs(c);
        if ac < ecp
            %entropy correction
            ac = (c^2+ecp^2)*0.5/ecp;
        end
        f = -(dt*c^2/dx*q+ac*(d0-q));
        flux(i) = 0.5*(c*ys4(i)+c*ys4(i+1)+f);
    end
    for i = 3:mx+1
        work(i) = ys4(i)-(dt/dx)*(flux(i)-flux(i-1));
    end
    for i = 3:mx+1
        ys4(i) = work(i);
    end
end
for i = 4:mx+2
    TV4 = TV4+abs(ys4(i)-ys4(i-1));
end

%% 1st order upwind scheme
TV5 = 0.0;
for i = 1:mx+3
    ys5(i) = yi(i);
end

for n = 1:nlast
    for i = 3:mx+1
        work(i) = ys5(i)-cfl*(ys5(i)-ys5(i-1));
    end
    for i = 3:mx+1
        ys5(i) = work(i);
    end
end
for i = 4:mx+2
    TV5 = TV5+abs(ys5(i)-ys5(i-1));
end

%% Lax-Wendroff scheme
TV6 = 0.0;
for i = 1:mx+3
    ys6(i) = yi(i);
end

for n = 1:nlast
    for i = 3:mx+1
        work(i) = 0.5*cfl*(1+cfl)*ys6(i-1)+(1-cfl^2)*ys6(i)-0.5*cfl*(1-cfl)*ys6(i+1);
    end
    for i = 3:mx+1
        ys6(i) = work(i);
    end
end
for i = 4:mx+2
    TV6 = TV6+abs(ys6(i)-ys6(i-1));
end
toc()
%% plot
figure(1);
plot(x(3:mx+2),yi(3:mx+2),x(3:mx+2),ye(3:mx+2),x(3:mx+2),ys1(3:mx+2),x(3:mx+2),ys2(3:mx+2),x(3:mx+2),ys3(3:mx+2),x(3:mx+2),ys4(3:mx+2));
ylim([-0.1 1.2]);
xlim([0.0 1.0]);
title("CFL="+string(cfl)+", after "+string(nlast)+"steps");
legend("Initial Condition","Strict","Minmod","Superbee","Sym-TVD(minmod)","Sym-TVD(superbee)")

figure(2);
plot(x(3:mx+2),yi(3:mx+2),x(3:mx+2),ye(3:mx+2),x(3:mx+2),ys3(3:mx+2),x(3:mx+2),ys4(3:mx+2),x(3:mx+2),ys5(3:mx+2),x(3:mx+2),ys6(3:mx+2));
ylim([-0.1 1.2]);
xlim([0.0 1.0]);
title("CFL="+string(cfl)+", after "+string(nlast)+"steps");
legend("Initial Condition","Strict","Sym-TVD(minmod)","Sym-TVD(superbee)","1st Upwind","Lax-Wendroff");

%% minmod
function valm = minmod(a,b)
valm = sign(a)*max(0.0,min(abs(a),b*sign(a)));
end


