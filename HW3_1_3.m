clc;clear;
f = @(t, x) t + 2*x*t;  % 함수 정의
x_exact = @(t) 0.5 * (exp(t.^2) - 1); % 실제 함수
t0 = 0; x0 = 0; % 초기 조건

t_end = 2;
h = 0.01;
N = (t_end - t0) / h; % 문제에 주어진 조건

t= t0:h:t_end;
x_rk2= zeros(1, N+1);  x_rk4= zeros(1, N+1);
x_rk2(1) = x0; x_rk4(1) = x0; %t,x 행렬 만들기

%2차 runge_kutta
for i = 1:N
    k1 = f(t(i), x_rk2(i));
    k2 = f(t(i) + h/2, x_rk2(i) + h*k1/2);
    x_rk2(i+1) = x_rk2(i) + h*k2;
end

%4차 runge_kutta
for j = 1:N
    k1 = f(t(j), x_rk4(j));
    k2 = f(t(j) + h/2, x_rk4(j) + h*k1/2);
    k3 = f(t(j) + h/2, x_rk4(j) + h*k2/2);
    k4 = f(t(j) + h, x_rk4(j) + h*k3);
    x_rk4(j+1) = x_rk4(j) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

%실제 함수
x_true=x_exact(t);

error2=max(abs(x_true-x_rk2))
error4=max(abs(x_true-x_rk4))
h2=error4/error2 %O(h^4)와 O(h^2)가 약 1/h^2배 만큼 차이나는 것을 확인

% 결과 그래프
plot(t, x_rk2, 'k-'); hold on;
plot(t, x_rk4, 'b-')
plot(t, x_true, 'r-')
xlabel('t');
ylabel('x(t)');
legend('x_r_k_2','x_r_k_4','x_t_r_u_e')
title('실제 함수와 rk2,rk4 비교');
grid on;