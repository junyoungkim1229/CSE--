clc;
clear;
f = @(t, x) t + 2*x*t;  % 함수 정의
t0 = 0; x0 = 0; % 초기 조건

t_end = 2;
h = 0.01;
N = (t_end - t0) / h; % 문제에 주어진 조건

t= t0:h:t_end;
x= zeros(1, N+1);
x(1) = x0; %t,x 행렬 만들기

% 2차 runge-kutta mid-point method
for n = 1:N
    k1 = f(t(n), x(n));
    k2 = f(t(n) + h/2, x(n) + h*k1/2);
    x(n+1) = x(n) + h*k2;
end

% 결과 그래프
plot(t, x, '-o');
xlabel('t');
ylabel('x(t)');
title('Second-order Runge-Kutta Method(midpoint)');
grid on;