%% Poisson: jacobi
clear; clc;

% 변수 설정
N     = 64; tol = 1e-8; maxit = 5e5;      
% 그리드
h = 1/(N+1); 
x = linspace(0,1,N+2);
y = linspace(0,1,N+2);

% f(x,y)
F = zeros(N+2,N+2);
for j = 1:N+2
    for i = 1:N+2
        F(i,j) = sin(pi*x(i))*sin(pi*y(j)); % 함수 정의
    end
end

% Initialization
u    = zeros(N+2,N+2);   % 변수 설정
unew = u;

% bnorm = ||f|| 
bnorm = 0.0;
for j = 2:N+1
    for i = 2:N+1
        bnorm = bnorm + F(i,j)^2;
    end
end
if bnorm == 0, bnorm = 1; end
bnorm = sqrt(bnorm);

res_hist = zeros(maxit,1); %잔차 행렬
t0 = tic; %계산 시간 측정 시작

% Jacobi iteration
for k = 1:maxit
    for j = 2:N+1
        for i = 2:N+1
            unew(i,j) = 0.25 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - (h^2)*F(i,j) );
        end   
    end

    % 잔차 계산
    res = 0.0;
    for j = 2:N+1
        for i = 2:N+1
            Lap_ij = (unew(i+1,j) + unew(i-1,j) + unew(i,j+1) + unew(i,j-1)- 4*unew(i,j) ) / (h^2);
            rij = F(i,j) - Lap_ij;
            res = res + rij^2;
        end
    end
    rrel = sqrt(res)/bnorm;
    res_hist(k) = rrel; %잔차 행렬
    u = unew;

    if rrel < tol
        res_hist = res_hist(1:k);
        break;
    end
end
t_elapsed = toc(t0); %계산 시간 측정 끝

% performance
fprintf('Jacobi (loops): N=%d, iters=%d, rel residual=%.3e, time=%.3fs\n', ...
        N, numel(res_hist), res_hist(end), t_elapsed);

%Plots (iteration,residual)
figure; semilogy(res_hist,'-o','LineWidth',1); grid on;
xlabel('iteration'); ylabel('relative residual');
title(sprintf('Jacobi residual history (N=%d)',N));

%Plots (u )
figure; 
surf(x,y,u'); shading interp; colorbar; view(30,30);
title('u(x,y)'); xlabel x; ylabel y;
