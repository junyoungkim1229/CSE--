clc;
clear;
n = 10;
i = 0:n;
x = cos((2*i + 1) * pi / (2*n + 2));  % Chebyshev node
y = 1 ./ (1 + 16 * x.^2);
X_int = linspace(-1, 1, 1000);
Ytrue = 1 ./ (1 + 16 * X_int.^2);
Yint = lagrange_interpolating(x, y, X_int);

plot(X_int, Ytrue, 'b', X_int, Yint, 'r--');
legend('f(x)', 'Lagrange p(x)');