clear; clc; close all;

% 주어진 조건
nu  = 1;         
n  = 2;        
U0 = 1;        
L  = 10;
Ny  = 101;                    
dy  = L/(Ny-1);                  
y   = linspace(0,L,Ny).';

t_final = 2*pi; % 비교 시점(한 주기 후)
lambda_max = 0.4;               

% Δt 세트 (FTCS 안정 조건 이내)
dt_base = lambda_max*dy^2/nu;    
dt_vec  = dt_base * (0.5).^(0:4);% 5 단계

% exact solution
eta  = @(yy) sqrt(n/(2*nu))*yy;
u_ex = @(yy,t) U0*exp(-eta(yy)).*cos(n*t - eta(yy));

% 결과 저장
err_FTCS = zeros(size(dt_vec));
err_CN   = zeros(size(dt_vec));

% 반복 : Δt 별 계산 
for kdt = 1:numel(dt_vec)
    dt = dt_vec(kdt);
    Nt = round(t_final/dt);  % 필요한 스텝 수
    lambda = nu*dt/dy^2;        
    
    % 1) FTCS
    u = zeros(Ny,1);
    for nt = 1:Nt
        t = (nt-1)*dt;
        u(1)   =  U0*cos(n*t);
        u(end) =  0;
        u(2:end-1) = u(2:end-1) + ...
            lambda*(u(3:end) - 2*u(2:end-1) + u(1:end-2));
    end
    err_FTCS(kdt) = sqrt(sum((u - u_ex(y,t_final)).^2)*dy);
    
    % 2) Crank–Nicolson
    Nint      = Ny-2;
    e         = ones(Nint,1);
    main_A    = (1+lambda)*e;
    off_A     = (-lambda/2)*e;
    A         = spdiags([off_A main_A off_A],[-1 0 1],Nint,Nint);
    [Lmat,Umat] = lu(A);                    % LU 1회
    B         = spdiags([-off_A (1-lambda)*e -off_A],[-1 0 1],Nint,Nint);
    
    u = zeros(Ny,1);                        % 초기화
    for nt = 1:Nt
        t = (nt-1)*dt;
        % 경계
        u(1)   =  U0*cos(n*t);
        u(end) =  0;
        % rhs
        rhs = B*u(2:end-1);
        rhs(1)   = rhs(1)   + (lambda/2)*(u(1)+u(1));
        rhs(end) = rhs(end) + (lambda/2)*(u(end)+u(end));
        u_int = Umat\(Lmat\rhs);
        u(2:end-1) = u_int;
    end
    err_CN(kdt) = sqrt(sum((u - u_ex(y,t_final)).^2)*dy);
end

% 그래프 출력
figure; loglog(dt_vec,err_FTCS,'-o','LineWidth',1.4); hold on;
loglog(dt_vec,err_CN  ,'-s','LineWidth',1.4);
grid on; xlabel('\Deltat'); ylabel('L_2 Error');
title('Time-step Convergence  (FTCS vs Crank–Nicolson)');
legend('FTCS (1^{st})','C–N (2^{nd})','Location','northwest');

% 기울기 계산
p_FTCS = polyfit(log(dt_vec),log(err_FTCS),1);
p_CN   = polyfit(log(dt_vec),log(err_CN),1);
fprintf('FTCS   수렴 차수 ≈ %.2f\n',-p_FTCS(1));
fprintf('C–N    수렴 차수 ≈ %.2f\n',-p_CN(1));
