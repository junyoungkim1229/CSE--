%% Linearity test for Poisson via SOR
clear; clc;

% 변수 설정
N     = 64;
tol   = 1e-5;
maxit = 5e5;
omega = 1.90;      

% 그리드
h = 1/(N+1);
x = linspace(0,1,N+2);
y = linspace(0,1,N+2);
[X,Y] = meshgrid(x,y);

% 함수 정의
F1 = sin(pi*X).*sin(pi*Y);
F2 = exp(-100*((X-0.5).^2 + (Y-0.5).^2));
F12 = F1 + F2;

% (1) f1+f2를 한 번에 계산한 u 계산
[u12, hist12] = sor_poisson(F12, N, h, tol, maxit, omega);

% (2) f2만으로 u2 계산
[u2,  hist2 ] = sor_poisson(F2,  N, h, tol, maxit, omega);

% (2)-1 f1만으로 u1 계산 (위에서 진행했던 것을 함수로만 정리)
[u1,  hist1 ] = sor_poisson(F1,  N, h, tol, maxit, omega);

% (3) 선형성 비교
lin_err_rel = norm(u12(2:end-1,2:end-1) - (u1(2:end-1,2:end-1)+u2(2:end-1,2:end-1)),'fro') / ...
              norm(u12(2:end-1,2:end-1),'fro');

fprintf('Linearity check:  ||u(f1+f2) - (u(f1)+u(f2))|| / ||u(f1+f2)|| = %.3e\n', lin_err_rel);

% 잔차 비교
figure; semilogy(hist1,'-','LineWidth',1); hold on;
semilogy(hist2,'--','LineWidth',1);
semilogy(hist12,'-.','LineWidth',1); grid on;
xlabel('iteration'); ylabel('relative residual');
legend('f1','f2','f1+f2'); title(sprintf('SOR(\\omega=%.2f) residual histories',omega));

% 각각 해 비교
figure; subplot(1,3,1); contourf(x,y,u1',20,'LineStyle','none'); colorbar; axis equal tight; title('u_1 (f_1)');
subplot(1,3,2); contourf(x,y,u2',20,'LineStyle','none'); colorbar; axis equal tight; title('u_2 (f_2)');
subplot(1,3,3); contourf(x,y,u12',20,'LineStyle','none'); colorbar; axis equal tight; title('u (f_1+f_2)');

figure; contourf(x,y,(u12-(u1+u2))',20,'LineStyle','none'); colorbar; axis equal tight;
title('u - (u_1+u_2)  (선형성 잔차)');

% SOR 함수(위에서 했던 코드를 함수로 만듦)
function [u, res_hist] = sor_poisson(F, N, h, tol, maxit, omega)
    u = zeros(N+2,N+2);  

    % ||f|| (interior)
    b2=0; 
    for j=2:N+1, 
        for i=2:N+1
            b2 = b2 + F(i,j)^2; 
        end
    end

    if b2==0, b2=1
    end 
    bnorm = sqrt(b2);

    res_hist = zeros(maxit,1);
    for k = 1:maxit
        for j = 2:N+1
            for i = 2:N+1
                uGS = 0.25*( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - h^2*F(i,j) );
                u(i,j) = (1-omega)*u(i,j) + omega*uGS;
            end
        end
        % residual r = f - Lap(u)
        r2=0;
        for j=2:N+1
            for i=2:N+1
                Lap = (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(h^2);
                rij = F(i,j) - Lap; r2 = r2 + rij^2;
            end
        end
        rrel = sqrt(r2)/bnorm; res_hist(k) = rrel;
        if rrel < tol
            res_hist = res_hist(1:k);
            break;
        end
    end
end
