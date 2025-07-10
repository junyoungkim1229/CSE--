function yspl = cubic_spline(x, y, x_fine)
    n = length(x);
    h = diff(x);
    

    A = zeros(n);    % n*n 행렬 생성  
    b = zeros(n, 1); % n*1행렬 생성

    % 공식 계수 입력
    for i = 2:n-1
        A(i,i-1) = h(i-1);
        A(i,i)   = 2 * (h(i-1) + h(i));
        A(i,i+1) = h(i);
        b(i) = 6 * ((y(i+1) - y(i)) / h(i) - (y(i) - y(i-1)) / h(i-1));
    end

    % 경계에서 2차 미분이 0이 되는 조건 적용
    A(1,1) = 1;
    A(n,n) = 1;
    b(1) = 0;
    b(n) = 0;
 
    % 이계도함수 찾기
    M = A \ b;

    % 스플라인 행렬의 크기 정의
    yspl = zeros(size(x_fine));

    for k = 1:length(x_fine)
        for i = 1:n-1
            if x_fine(k) >= x(i) && x_fine(k) <= x(i+1) %x_fine이 속하는 i값 찾기
                hi = x(i+1) - x(i);
                a = (x(i+1) - x_fine(k)) / hi;
                b = (x_fine(k) - x(i)) / hi;
                yspl(k) = a*y(i) + b*y(i+1) + ((a^3 - a)*M(i) + (b^3 - b)*M(i+1)) * (hi^2) / 6;
                break;
            end
        end
    end
end