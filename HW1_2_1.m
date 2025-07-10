clc;
clear;
f = @(x) x.^5 - 9*x.^4 - x.^3 + 17*x.^2 - 8*x - 8; % 함수 정의
[A1,A2]=bisection_method(f,-10,-1,10^-8);          % -10 ~ -1  구간해
[B1,B2]=bisection_method(f,-1,0,10^-8);            % -1  ~  0  구간해
[C1,C2]=bisection_method(f,0,10,10^-8);            %  0  ~  10 구간해
ans=[A1 A2;B1 B2;C1 C2]

