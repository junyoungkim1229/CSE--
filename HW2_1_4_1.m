clc;
clear;
f = @(x) 1 ./ (1 + 16*x.^2); %함수 입력
x = -1:0.2:1;                % uniform nodes
y = f(x);
x_fine = linspace(-1, 1, 1000);  %원래 함수 좌표

% Natural cubic spline 보간
y_spline = cubic_spline(x, y, x_fine);

% Plotting
plot(x_fine, f(x_fine), 'b', 'LineWidth', 1.5); hold on;
plot(x_fine, y_spline, 'r--', 'LineWidth', 1.5);
legend('f(x)', 'Cubic Spline');
xlabel('x'); ylabel('y');
grid on;