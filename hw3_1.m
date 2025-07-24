clear; clc; close all;

% 문제에 주어진 조건
nu  = 1;                
n   = 2;                
U0  = 1;                
L   = 10;              
Ny  = 201;              % 점의 개수 
dy  = L/(Ny-1);         % 간격
y   = linspace(0,L,Ny).';   % y 벡터 (열)

lambda = 0.4;           % ν·dt/dy²  (안정 조건 λ≤0.5)
dt = lambda*dy^2/nu;    % time step
fprintf('λ = %.4f (OK ≤ 0.5)\n',lambda);

% 시간 설정
T_transient = 10*pi;                   % transient 시간
t_phase     = [0, pi/2, pi, 3*pi/2, 2*pi];   % transient 구간 시점
t_phase_qs  = T_transient + t_phase;     % quasi-steady 시점
t_max       = t_phase_qs(end);           % 최종 시간
Nt          = round(t_max/dt)+1;         % 총 스텝 수
t_vec       = (0:Nt-1)*dt;               % 시간 배열

idx_tr   = round(t_phase   /dt)+1;       % 인덱스(transient)
idx_qs   = round(t_phase_qs/dt)+1;       % 인덱스(quasi-steady)

% 행렬 생성
u      = zeros(Ny,1);
u_new  = zeros(Ny,1);
utr    = cell(numel(t_phase),1);
uqs    = cell(numel(t_phase),1);

% FTCS 구현
for k = 1:Nt
    t = t_vec(k);

    % 경계 조건
    u(1)   =  U0*cos(n*t);   % y=0
    u(end) = 0;              % y=L

    % 내부 노드 업데이트 (FTCS)
    for i = 2:Ny-1
        u_new(i) = u(i) + lambda*(u(i+1) - 2*u(i) + u(i-1));
    end
    u = u_new;   % 다음 스텝 준비

    % 시점 저장
    idx_now = k;
    tr_hit  = find(idx_tr  == idx_now,1);
    qs_hit  = find(idx_qs  == idx_now,1);
    if ~isempty(tr_hit), utr{tr_hit} = u; 
    end
    if ~isempty(qs_hit), uqs{qs_hit} = u; 
    end
end

% Transient 그래프
figure('Name','Transient'); hold on;
clr = lines(numel(t_phase));
for i = 1:numel(t_phase)
    plot(utr{i}, y, 'Color', clr(i,:), ...
        'DisplayName', sprintf('t = %.2f', t_phase(i)));
end
xlabel('u'); ylabel('y'); grid on;
title('Transient Velocity Profiles');
legend show;

% Quasi-Steady 그래프
figure('Name','Quasi-Steady'); hold on;
for i = 1:numel(t_phase)
    plot(uqs{i}, y, 'Color', clr(i,:), ...
        'DisplayName', sprintf('t = %.2f', t_phase(i)));
end
xlabel('u'); ylabel('y'); grid on;
title('Quasi-Steady Velocity Profiles  (T = 10π 이후)');
legend show;
