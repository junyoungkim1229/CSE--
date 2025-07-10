clc;
clear;
f = @(x) 1 ./ (1 + 16*x.^2); %함수 입력
n = 10;
i=0:n; % Chebyshev 노드
x_cheb = cos((2*i + 1) * pi / (2*n + 2));  % [-1, 1] 구간
x_cheb=sort(x_cheb); % 순서 정렬
y_cheb = 1 ./ (1 + 16 * x_cheb.^2);

x_fine = linspace(-1, 1, 1000); %원래 함수 좌표

% 자연 삼차 스플라인 보간 수행
y_cheb_spline = cubic_spline(x_cheb, y_cheb, x_fine);

% 그래프 출력
plot(x_fine, f(x_fine), 'b', 'LineWidth', 1.5); hold on;
plot(x_fine, y_cheb_spline, 'r--', 'LineWidth', 1.5);
legend('f(x)', 'Cubic Spline');
xlabel('x'); ylabel('y');
grid on;