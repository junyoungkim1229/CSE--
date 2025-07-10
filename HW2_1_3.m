clc;
clear
x_plot = linspace(-1, 1, 1000);

% uniform 노드
x_uniform = -1:0.2:1;
w_uniform = ones(size(x_plot));
k=length(x_uniform);
for i = 1:k
    w_uniform = w_uniform .* abs(x_plot - x_uniform(i));
end

% Chebyshev 노드 (1st kind)
n = 10;
i = 0:n;
x_cheb = cos((2*i + 1) * pi / (2*(n+1)));
w_cheb = ones(size(x_plot));
m=length(x_cheb);
for i = 1:m
    w_cheb = w_cheb .* abs(x_plot - x_cheb(i));
end

% 그래프 비교
plot(x_plot, w_uniform, 'r--', 'LineWidth', 1.5); hold on;
plot(x_plot, w_cheb, 'b-', 'LineWidth', 1.5);
legend('Uniform nodes', 'Chebyshev nodes');
xlabel('x'); ylabel('∏|x - x_i|');
grid on;