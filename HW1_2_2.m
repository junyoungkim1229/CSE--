clc;
clear;
f = @(x) x.^5 - 9*x.^4 - x.^3 + 17*x.^2 - 8*x - 8; % 함수 정의
df= @(x) 5*x.^4 - 36*x.^3 - 3*x.^2 + 34*x - 8;     % 도함수 정의
[A1,A2]=newton_method(f,df,-10,10^-8,10000);       % x0=-10의 해
[B1,B2]=newton_method(f,df,-0.1,10^-8,10000);      % x0=-0.1의 해
[C1,C2]=newton_method(f,df,10,10^-8,10000);        % x0=10의 해
ans=[A1 A2;B1 B2;C1 C2]
