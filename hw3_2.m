clear; clc; close all;

% 문제에 주어진 조건
nu  = 1;     
n  = 2;     
U0 = 1;     
L  = 10;
Ny  = 201;                     
dy  = L/(Ny-1);                
y   = linspace(0,L,Ny).';     
dt  = 0.01;                   
lambda = nu*dt/dy^2;           

% 시간 설정 
T_tr = 10*pi;                          % transient 기간
tP   = [0, pi/2, pi, 3*pi/2, 2*pi];    % 5개 관측 시점
tP_q = T_tr + tP;                      % quasi-steady 시점
t_end = tP_q(end);
Nt    = round(t_end/dt)+1;
t_vec = (0:Nt-1)*dt;

idx_tr = round(tP   /dt)+1;   % 인덱스(transient)
idx_q  = round(tP_q /dt)+1;   % 인덱스(quasi-steady)

% C-N 계수행렬
Nint = Ny-2;                     % 내부 노드 수

% Crank–Nicolson 계수행렬 A (좌변)
e     = ones(Nint,1);            % 편의용 단위벡터
main  = (1+lambda)*e;            % 길이 Nint
off   = (-lambda/2)*e;           % 길이 Nint

A = spdiags([off main off],[-1 0 1],Nint,Nint);   % 좌변 행렬
[Lmat,Umat] = lu(A);                              % LU 분해
B = spdiags([ -off (1-lambda)*e -off ],[-1 0 1],Nint,Nint);  % 우변 행렬



% 변수 
u     = zeros(Ny,1);     % 이전 스텝
u_new = u;               % 다음 스텝
utr  = cell(numel(tP),1);
uqs  = cell(numel(tP),1);

% 시간 적분
for k = 1:Nt
    t = t_vec(k);

    % 경계조건
    u(1)   =  U0*cos(n*t);   % y=0
    u(end) = 0;              % y=L

    % 우변 벡터 계산 (내부 노드만)
    rhs = B*u(2:end-1);

    % 경계 조건이 rhs에 주는 영향
    rhs(1)   = rhs(1)   + (lambda/2)*(u(1)      + u(1));        % 좌측
    rhs(end) = rhs(end) + (lambda/2)*(u(end)    + u(end));      % 우측

    % 선형계 풀기 
    ytmp = Lmat\rhs;
    u_int_new = Umat\ytmp;

    % 결과
    u_new(2:end-1) = u_int_new;
    u = u_new;

    % 저장
    hit_tr = find(idx_tr==k,1);  if ~isempty(hit_tr), utr{hit_tr}=u; end
    hit_q  = find(idx_q ==k,1);  if ~isempty(hit_q ), uqs{hit_q }=u; end
end

% 그래프 출력
clr = lines(numel(tP));

figure('Name','Transient (C-N)'); hold on;
for i = 1:numel(tP)
    plot(utr{i}, y,'Color',clr(i,:), ...
        'DisplayName',sprintf('t = %.2f',tP(i)));
end
xlabel('u'); ylabel('y'); grid on; title('Transient Velocity (C-N)');
legend show;

figure('Name','Quasi-Steady (C-N)'); hold on;
for i = 1:numel(tP)
    plot(uqs{i}, y,'Color',clr(i,:), ...
        'DisplayName',sprintf('t = %.2f',tP(i)));
end
xlabel('u'); ylabel('y'); grid on; title('Quasi-Steady Velocity (C-N)');
legend show;
