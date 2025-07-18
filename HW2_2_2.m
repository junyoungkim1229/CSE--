clc; clear;

f = @(x) sin(16 - x.^2);  % 함수 정의
f_exact = @(x) -2 * cos(16 - x.^2) - 4 * x.^2 .* sin(16 - x.^2);  %이계도함수 정의
h = 8 / 32;  % 33개 점을 사용하므로 간격을 8/(n-1)로 설정

x = 0:h:8; %점 생성
x_fine=linspace(0,8,1000); %원래 함수 표현을 위한 점 생성

% second-order central difference scheme
f2_approx = zeros(1, length(x)); %근사할 이계도함수 행렬 생성
for i = 2:length(x)-1
    f2_approx(i) = (f(x(i+1)) - 2*f(x(i)) + f(x(i-1))) / h^2;
end

% 정확한 두 번째 도함수 계산
f2_exact = f_exact(x_fine);

% 그래프 그리기
plot(x, f2_approx, 'r-o'); hold on;
plot(x_fine, f2_exact, 'b-');
legend('Approximate f''''(x)','Exact f''''(x)');
xlabel('x'); ylabel('f''''(x)');
title('비교');
grid on;