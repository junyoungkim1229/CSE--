clear; clc; close all;

% 문제에 주어진 조건
nu = 1;  n  = 2;  U0 = 1;  L = 2;          % L=2로 변경
Ny = 161;                
dy = L/(Ny-1);    
y  = linspace(0,L,Ny).';
dt = 0.006;               
lambda = nu*dt/dy^2;   % C–N이므로 안정
fprintf('lambda = %.3f\n',lambda);

% 시간 설정
T_tr = 10*pi;                         % transient 기간
tP   = [0, pi/2, pi, 3*pi/2, 2*pi];   % 5개 시점
tP_q = T_tr + tP;   t_end = tP_q(end);
Nt   = round(t_end/dt)+1;  t_vec=(0:Nt-1)*dt;
idxP = round(tP  /dt)+1;   idxQ = round(tP_q/dt)+1;

% Crank–Nicolson 계수행렬
Nint = Ny-2;  e = ones(Nint,1);
main = (1+lambda)*e;      off = (-lambda/2)*e;
A = spdiags([off main off],[-1 0 1],Nint,Nint);
[Lm,Um] = lu(A);
B = spdiags([-off (1-lambda)*e -off],[-1 0 1],Nint,Nint);

% 적분
u=zeros(Ny,1);  u_new=u;
utr=cell(numel(tP),1); uqs=cell(numel(tP),1);

for k=1:Nt
    t = t_vec(k);
    u(1)=U0*cos(n*t); u(end)=0;                     % 경계 조건
    rhs = B*u(2:end-1);
    rhs(1)   = rhs(1)+ (lambda/2)*(u(1)+u(1));
    rhs(end) = rhs(end)+(lambda/2)*(u(end)+u(end));
    u(2:end-1)= Um\(Lm\rhs);                        % 해
    % 저장
    hit=find(idxP==k); if ~isempty(hit),  utr{hit}=u; end
    hit=find(idxQ==k); if ~isempty(hit),  uqs{hit}=u; end
end

% 그래프 작성
clr = lines(numel(tP));
figure('Name','Transient L=2'); hold on;
for i=1:numel(tP), plot(utr{i},y,'Color',clr(i,:), ...
      'DisplayName',sprintf('t=%.2f',tP(i))); end
xlabel('u'); ylabel('y'); title('Transient  (L=2)');
legend show; grid on;

figure('Name','Quasi-Steady L=2'); hold on;
for i=1:numel(tP), plot(uqs{i},y,'Color',clr(i,:), ...
      'DisplayName',sprintf('t=%.2f',tP(i))); end
xlabel('u'); ylabel('y'); title('Quasi-Steady  (L=2,  T=10\pi 이후)');
legend show; grid on;
