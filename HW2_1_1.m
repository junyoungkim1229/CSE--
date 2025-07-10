clc;
clear;
x = -1:0.2:1;
y = 1 ./ (1 + 16 * x.^2);
X_int = linspace(-1, 1, 1000);
Ytrue = 1 ./ (1 + 16 * X_int.^2);
Yint = lagrange_interpolating(x, y, X_int);

plot(X_int, Ytrue, 'b', X_int, Yint, 'r--');
legend('f(x)', 'Lagrange p(x)');