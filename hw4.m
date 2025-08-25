clear; clc; close all;

%% -------- 기본 설정 --------
a        = 1.0;             % thermal diffusivity
tol_ss   = 1e-10;            % 정상상태(고정점) 판정 기준: ||u^{n+1}-u^n||_inf
maxSteps = 2e5;              % 안전장치
Mser     = 61; Nser = 61;    % 정해(정상상태) 급수의 모드 수 (합리적: 31~101, 홀/짝은 내부에서 처리)

% (A) 기본 케이스(표면 비교/오차 시각화)
Nx0 =61; Ny0 = 61; dt0 = 5e-2;

% (B) Δt 변화에 따른 수렴속도/오차
dt_list = [2e-3 5e-3 1e-2 2e-2 5e-2];

% (C) 격자(공간) 수렴 (Ny=Nx로 동일 증가)
Nx_list = [16 32 48 64 96];  Ny_list = Nx_list;
dt_for_hconv = 1e-2;

%% -------- (A) 수치해 vs 정해(정상상태) --------
[xi, yi, X, Y, Sxy] = make_grid_and_source(Nx0,Ny0);

[U_num, steps] = crank_nicolson_to_steady(Nx0,Ny0,dt0,a,Sxy,tol_ss,maxSteps);
U_num = reshape(U_num, Nx0, Ny0).';      % (Ny x Nx)

U_ex  = exact_steady_on_grid(xi, yi, Mser, Nser); % (Ny x Nx)
E     = U_num - U_ex;

fprintf('[Base] Nx=%d, Ny=%d, dt=%.3g -> steps=%d,  ||err||_inf=%.3e\n', ...
    Nx0,Ny0,dt0,steps,norm(E(:),inf));

figure(1); 
surf(X,Y,U_num,'EdgeColor','none'); view(35,30); colorbar;
xlabel('x'); ylabel('y'); zlabel('\phi_{num}');
title(sprintf('Numerical steady state (CN), N=(%d,%d), dt=%.3g',Nx0,Ny0,dt0));

figure(2);
surf(X,Y,U_ex,'EdgeColor','none'); view(35,30); colorbar;
xlabel('x'); ylabel('y'); zlabel('\phi_{exact,ss}');
title(sprintf('Exact steady-state series, modes (%d,%d)',Mser,Nser));

figure(3);
surf(X,Y,abs(E),'EdgeColor','none'); view(35,30); colorbar;
xlabel('x'); ylabel('y'); zlabel('|error|');
title('Absolute error: | \phi_{num} - \phi_{exact,ss} |');

%% -------- (B) Δt 변화: 수렴속도(스텝 수) & 오차 --------
steps_dt = zeros(size(dt_list));
err_dt   = zeros(size(dt_list));

U_ex_base = U_ex;  % 동일 격자에서의 정해-정상상태

for k = 1:numel(dt_list)
    dt = dt_list(k);

    % 수치해 계산
    [Uvec, steps] = crank_nicolson_to_steady(Nx0,Ny0,dt,a,Sxy,tol_ss,maxSteps);
    U = reshape(Uvec, Nx0, Ny0).';  % (Ny x Nx)

    % 오차
    tmp = U - U_ex_base;
    err_dt(k)   = norm(tmp(:), inf);
    steps_dt(k) = steps;

    fprintf('dt=%.3g: steps=%d, ||err||_inf=%.3e\n', dt, steps, err_dt(k));
end



figure(4);
loglog(dt_list, max(err_dt,eps), 'o-','LineWidth',1.5); grid on;
xlabel('\Delta t'); ylabel('||error||_\infty');
title('Error vs \Delta t (steady state)');

%% ---- 시간 차수 검증 (고정 t*에서) ----
NxT = 96; NyT = 96;          % 충분히 촘촘하게
[xiT, yiT, ~, ~, SxyT] = make_grid_and_source(NxT,NyT);
tstar = 0.2;                  % 비교할 고정 시간
dt_listT = [1e-3 2e-3 5e-3 1e-2 2e-2];  % 여러 Δt

% t*의 정해(1번에서 유도한 시간 의존 해) - 급수 모드 수
Mser = 61; Nser = 61;
U_exact_t = exact_time_on_grid(xiT, yiT, tstar, Mser, Nser);  % (NyT x NxT)

err_dt = zeros(size(dt_listT));
for k=1:numel(dt_listT)
    dt = dt_listT(k);
    Uvec = crank_nicolson_to_time(NxT,NyT,dt,1.0,SxyT,tstar); % 아래 함수 추가
    U = reshape(Uvec, NxT, NyT).';
    tmp = U - U_exact_t;
    err_dt(k) = norm(tmp(:), inf);
    fprintf('dt=%.3g  ->  ||err(t*)||_inf = %.3e\n', dt, err_dt(k));
end

figure(6); 
loglog(dt_listT, max(err_dt,eps),'o-','LineWidth',1.5); grid on;
xlabel('\Delta t'); ylabel('||error at t^*||_\infty');
title('Temporal order at fixed t^*');

p_time = polyfit(log(dt_listT), log(max(err_dt,eps)), 1);
fprintf('Observed temporal order ~ %.3f (expected ~2 for CN)\n', -p_time(1));


%% -------- (C) 공간 수렴: ||err|| vs h --------
err_h = zeros(size(Nx_list));
hlist = zeros(size(Nx_list));

for k = 1:numel(Nx_list)
    Nx = Nx_list(k); Ny = Ny_list(k);
    [xi, yi, ~, ~, Sxy] = make_grid_and_source(Nx,Ny);
    [U, steps] = crank_nicolson_to_steady(Nx,Ny,dt_for_hconv,a,Sxy,tol_ss,maxSteps); %#ok<ASGLU>
    U = reshape(U, Nx, Ny).';
    Uex = exact_steady_on_grid(xi, yi, Mser, Nser);

    hx = 2/(Nx+1); hy = 2/(Ny+1);
    hlist(k) = max(hx,hy);
    temp1=U - Uex;
    err_h(k) = norm(temp1(:), inf);
    fprintf('N=(%d,%d): h=%.4f, ||err||_inf=%.3e\n', Nx,Ny,hlist(k),err_h(k));
end

figure(5);
loglog(hlist, max(err_h,eps), 'o-','LineWidth',1.5); grid on;
xlabel('h (max(h_x,h_y))'); ylabel('||error||_\infty');
title('Spatial convergence (expect ~O(h^2))');

%% ===================== 보조 함수들 =====================

function [xi, yi, X, Y, Sxy] = make_grid_and_source(Nx,Ny)
    % 내부 격자점 (Dirichlet 0 경계는 제외)
    hx = 2/(Nx+1);  hy = 2/(Ny+1);
    xi = (-1+hx) : hx : (1-hx);               % 1 x Nx
    yi = (-1+hy) : hy : (1-hy);               % 1 x Ny
    [X,Y] = meshgrid(xi, yi);                 % (Ny x Nx)
    Sxy = 2*(2 - X.^2 - Y.^2);                % source on interior grid
end

function [u, steps] = crank_nicolson_to_steady(Nx,Ny,dt,a,Sxy,tol_ss,maxSteps)
    % Crank–Nicolson: (I - dt/2*A) u^{n+1} = (I + dt/2*A) u^n + dt*b
    % A = a*(kron(Iy,Lx) + kron(Ly,Ix)), Dirichlet 0 (interior unknowns만)
    hx = 2/(Nx+1);  hy = 2/(Ny+1);
    ex = ones(Nx,1); ey = ones(Ny,1);

    Lx = spdiags([ex -2*ex ex], -1:1, Nx, Nx) / hx^2;
    Ly = spdiags([ey -2*ey ey], -1:1, Ny, Ny) / hy^2;

    Ix = speye(Nx); Iy = speye(Ny);
    A  = a*(kron(Iy,Lx) + kron(Ly,Ix));   % size (Nx*Ny) x (Nx*Ny)

    Ntot = Nx*Ny; I = speye(Ntot);
    LHS = (I - 0.5*dt*A);
    RHS = (I + 0.5*dt*A);

    b = Sxy(:);                    % time-independent source
    u = zeros(Ntot,1);
    steps = 0;

    % 미리 분해(반복 루프에서 빠르게)
    [Lfac, Ufac, P, Q] = lu(LHS);  % sparse LU

    for n = 1:maxSteps
        rhs = RHS*u + dt*b;
        u_new = Q*(Ufac \ (Lfac \ (P*rhs)));
        steps = n;
        if norm(u_new - u, inf) < tol_ss
            u = u_new; break;
        end
        u = u_new;
    end
end

function Uex = exact_steady_on_grid(xi, yi, Mser, Nser)
    % 정상상태 정해 (t->∞):
    % phi_ss = sum_{m,n>=1} [ S_mn / ( (pi^2/4)*(m^2+n^2) ) ] * sin(mπ(x+1)/2)*sin(nπ(y+1)/2)
    % where S_mn = (64/π^4)*(1-(-1)^m)*(1-(-1)^n)*(1/(m^3 n)+1/(m n^3))
    [X,Y] = meshgrid(xi, yi);
    Uex = zeros(size(X));
    pi2 = pi^2;  % just for speed
    for m = 1:Mser
        fm = 1 - (-1)^m;   % 0 (짝수) or 2 (홀수)
        if fm==0, continue; end
        sx = sin(m*pi*(X+1)/2);
        for n = 1:Nser
            fn = 1 - (-1)^n;
            if fn==0, continue; end
            Smn = (64/pi^4) * (fm*fn) * (1/(m^3*n) + 1/(m*n^3));
            coeff = Smn / ( (pi2/4)*(m^2+n^2) );
            Uex = Uex + coeff * sx .* sin(n*pi*(Y+1)/2);
        end
    end
end

function U = exact_time_on_grid(xi, yi, t, Mser, Nser)
    [X,Y] = meshgrid(xi, yi);
    U = zeros(size(X));
    for m=1:Mser
        fm = 1-(-1)^m; if fm==0, continue; end
        sx = sin(m*pi*(X+1)/2);
        for n=1:Nser
            fn = 1-(-1)^n; if fn==0, continue; end
            Smn  = (64/pi^4)*(fm*fn)*(1/(m^3*n)+1/(m*n^3));
            lam  = (pi^2/4)*(m^2+n^2);
            coeff= Smn/lam * (1 - exp(-lam*t));
            U = U + coeff * sx .* sin(n*pi*(Y+1)/2);
        end
    end
end

function u = crank_nicolson_to_time(Nx,Ny,dt,a,Sxy,tstar)
    hx = 2/(Nx+1); hy = 2/(Ny+1);
    ex = ones(Nx,1); ey = ones(Ny,1);
    Lx = spdiags([ex -2*ex ex], -1:1, Nx, Nx) / hx^2;
    Ly = spdiags([ey -2*ey ey], -1:1, Ny, Ny) / hy^2;
    Ix = speye(Nx); Iy = speye(Ny);
    A  = a*(kron(Iy,Lx) + kron(Ly,Ix));
    Ntot = Nx*Ny; I = speye(Ntot);
    LHS = (I - 0.5*dt*A); RHS = (I + 0.5*dt*A);
    b = Sxy(:);
    u = zeros(Ntot,1);

    nsteps = ceil(tstar/dt);
    [Lfac,Ufac,P,Q] = lu(LHS);
    for n=1:nsteps
        rhs = RHS*u + dt*b;
        u = Q*(Ufac\(Lfac\(P*rhs)));
    end
end
