clear
close all
% 基本參數
a_p = -5;
k_p = 3 ;
u   = 1 ;
% 設定X軸畫圖間隔大小，以每0.001秒畫一次
h   = 0.001;
% 總共取3秒為畫圖區間，間隔為h
t   = 0:h:3; 
% 設定微分方程初始值x(t=0)=0
x1   = 0;
% 利用backward difference計算微分，x'=(x(t)-x(t-h))/h
% 原式可改寫成(x(t)-x(t-h))/h = ax(t)+ku，移項整理後使用迴圈將x(t)求出
for ii = 2:length(t)
    x1(ii) = (x1(ii-1) + h*k_p*u)/(1-a_p*h);
end
% 將x對t畫圖
figure
plot(t,x1)
% 依序設置圖片標題、命名X軸、命名Y軸
title('u=1')
xlabel('time')
ylabel('x_p')
% 變更u函數
u  = sin(t);
x2 = 0;
for ii = 2:length(t)
    x2(ii) = (x2(ii-1) + h*k_p*u(ii))/(1-a_p*h);
end
figure
plot(t,x2)
title('u=sint')
xlabel('time')
ylabel('x_p')
% 變更u函數
u  = t.^2;
x3 = 0;
for ii = 2:length(t)
    x3(ii) = (x3(ii-1) + h*k_p*u(ii))/(1-a_p*h);
end
figure
plot(t,x3)
title('u=t^2')
xlabel('time')
ylabel('x_p')
% 變更u函數
u  = t.^2.*sin(t);
x4 = 0;
for ii = 2:length(t)
    x4(ii) = (x4(ii-1) + h*k_p*u(ii))/(1-a_p*h);
end
figure
plot(t,x4)
title('u=t^2sint')
xlabel('time')
ylabel('x_p')
% 將不同的u繪製在同一張圖
figure
plot(t,x1,t,x2,t,x3,t,x4)
legend('u=1','u=sint','u=t^2','u=t^2sint')
xlabel('time')
ylabel('x_p')
