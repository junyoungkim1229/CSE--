clc;
clear;
f = @(x) sin(16 - x.^2);  % 함수 정의
f_exact = @(x) -2 * cos(16 - x.^2) - 4 * x.^2 .* sin(16 - x.^2);  % 이계도함수

n = [33, 65, 129]; %주어진 n값 행렬로 생성
errors = [];

for N = n
    h = 8 / (N - 1);  % 격자 간격
    x = linspace(0, 8, N);  % 격자점
    f2_approx = zeros(1, N);
    for i = 2:N-1
        f2_approx(i) = (f(x(i+1)) - 2*f(x(i)) + f(x(i-1))) / h^2; %second-order central difference scheme
    end

    f2_exact_values = f_exact(x); %실제 이계도함수 함숫값

    error = sqrt(sum((f2_exact_values(2:N-1) - f2_approx(2:N-1)).^2)*h);  % 첫번째 인덱스와 마지막 인덱스를 제외한 L2 norm error
    errors = [errors, error];
end

% x,y축 값을 각각 로그 스케일로 변경
log_dx = log10(8 ./ (n - 1));  % x축: log(Δx)
log_error = log10(errors);  % y축: log(error)

% 그래프 그리기
plot(log_dx, log_error, 'o-');
xlabel('log(\Delta x)');
ylabel('log(\epsilon)');
title('오차분석');
grid on;