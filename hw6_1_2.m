%% Poisson: gauss-seidel
clear; clc;

% 변수 설정
N     = 64;        % interior points per direction
tol   = 1e-8;      % relative residual tolerance
maxit = 5e5;       % max iterations

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
u = zeros(N+2,N+2);   % in-place update

% bnorm = ||f|| (interior)
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

% Gauss-Seidel iteration (좌->우, 하->상 스윕)
for k = 1:maxit
    for j = 2:N+1
        for i = 2:N+1
            u(i,j) = 0.25 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - (h^2)*F(i,j) );
        end %갱신된 값을 바로 사용
    end

    % 잔차 = f - Lap(u)
    res = 0.0;
    for j = 2:N+1
        for i = 2:N+1
            Lap_ij = (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(h^2);
            rij = F(i,j) - Lap_ij;
            res = res + rij^2;
        end
    end
    rrel = sqrt(res)/bnorm;
    res_hist(k) = rrel; %잔차 행렬

    if rrel < tol
        res_hist = res_hist(1:k);
        break;
    end
end
t_elapsed = toc(t0); %계산 시간 측정 끝

% performance
fprintf('gauss–Seidel method: N=%d, iters=%d, rel residual=%.3e, time=%.3fs\n', ...
        N, numel(res_hist), res_hist(end), t_elapsed);

% Plots
figure; semilogy(res_hist,'-o','LineWidth',1); grid on;
xlabel('iteration'); ylabel('relative residual');
title(sprintf('Gauss-Seidel residual history (N=%d)',N));

figure;
surf(x,y,u'); shading interp; colorbar; view(30,30);
title('u(x,y) - Gauss-Seidel'); xlabel x; ylabel y;
