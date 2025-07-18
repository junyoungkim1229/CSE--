clc; clear;format short

f = @(t, x) t + 2*x*t;  % 함수 정의
x_exact = @(t) 0.5 * (exp(t.^2) - 1); % 실제 함수
t0 = 0; x0 = 0; % 초기 조건

t_end = 2;
h_value= [0.01, 0.05, 0.1];
error = zeros(1, length(h_value));

for j = 1:3
    h = h_value(j);
    N = (t_end - t0) / h; % 문제에 주어진 조건
    t = t0:h:t_end;
    x = zeros(1, N+1);
    x(1) = x0; %t,x 행렬 만들기

    %4차 runge-kutta 
    for n = 1:N
    k1 = f(t(n), x(n));
    k2 = f(t(n) + h/2, x(n) + h*k1/2);
    k3 = f(t(n) + h/2, x(n) + h*k2/2);
    k4 = f(t(n) + h, x(n) + h*k3);
    x(n+1) = x(n) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end


%실제 함수
x_true=x_exact(t);

error(j)=abs(x_true(end)-x(end));
end

error_ratio=[error(2)/error(1) error(3)/error(2)] 
%각각 h가 5배, 2배 되었으므로, O(h^4)는 625배,16배 되어야한다.

figure;
plot(h_value, error, 'o-');
xlabel('Step size (h)');
ylabel('Error');
title('Error vs Step size (h) for 4th Order Runge-Kutta');
grid on;