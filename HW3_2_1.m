clc; clear;

%문제에 주어진 조건 설정
g = 9.81;                % 중력 가속도 (m/s^2)
B = 4.1e-4;              % Magnus 힘 계수
omega = 1800 * 2*pi / 60;  % 회전속도 rad/s (1800 rpm을 rad/s로 변환)
phi = deg2rad(135);       % 회전 방향 각도 (phi)

% 초기 속도 및 조건
v0 = 30;                 % 초기 속도 (m/s)
theta = deg2rad(1);      % 발사 각도
vx0 = v0 * cos(theta);   % x축 속도
vy0 = 0;                 % y축 속도
vz0 = v0 * sin(theta);   % z축 속도
x0 = 0; y0 = 0; z0 = 1.7; % 초기 위치 (지상에서 1.8m)

% 시간 설정
h = 0.001;               % 시간 간격
t_max = 5;               % 최대 시간
t = 0:h:t_max;

% 초기 조건 행렬
state = zeros(6, length(t));
state(:,1) = [x0; y0; z0; vx0; vy0; vz0];

% F(V) 함수 정의 (drag)
F = @(V) 0.0039 + 0.0058 ./ (1 + exp((V - 35) / 5));

% 4차 Runge-Kutta
for n = 1:length(t)-1
    s = state(:,n);
    k1 = baseball_rhs(s, F, B, g, phi, omega);
    k2 = baseball_rhs(s + 0.5 * k1, F, B, g, phi, omega);
    k3 = baseball_rhs(s + 0.5 * k2, F, B, g, phi, omega);
    k4 = baseball_rhs(s + k3, F, B, g, phi, omega);
    state(:,n+1) = s + (k1 + 2*k2 + 2*k3 + k4)*h/6;
    
    % if문을 통해 x > 18.39m 도달 시 정지하도록 설정
    if state(1,n+1) >= 18.39
        state = state(:,1:n+1); % 잘라냄
        t = t(1:n+1);
        break;
    end
end

% 결과 시각화
figure;
plot3(state(1,:), state(2,:), state(3,:), 'LineWidth', 2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Baseball Trajectory');
grid on;
view(3);


