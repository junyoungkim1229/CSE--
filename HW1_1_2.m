clc;
clear;
format short
A=randi([1,10],4,4);
B=randi([1,10],4,4);
C=A*B;%행렬곱

D=[C eye(4)] %eye: 단위행렬 생성
for i = 1:4    
    if D(i,i) == 0
        for j = i+1:4
            if D(j,i) ~= 0
                F = D(i,:);
                D(i,:) = D(j,:);
                D(j,:) = F; %자리 바꾸기
                break; %j 반복문 정지
            end
        end
    end %사다리꼴 행렬 만들기
    D(i,:) = D(i,:) / D(i,i); %피벗 계수 1로 만들기
    for j = 1:4
        if j ~= i
            D(j,:) = D(j,:) - D(j,i) * D(i,:);%가우스 소거법
        end
    end
end
C_inv=D(:,5:8)
I=C*C_inv