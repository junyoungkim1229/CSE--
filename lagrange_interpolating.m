function Yint = lagrange_interpolating(x, y, Xint)
    n = length(x);
    m = length(Xint);
    Yint = zeros(size(Xint));
    for k=1:m
        for i = 1:n
            L(i) = 1;
            for j = 1:n
                if j ~= i
                    L(i) = L(i) * (Xint(k) - x(j)) / (x(i) - x(j));
                end
            end  
        end  
        Yint(k) = sum(y .* L);
    end
  