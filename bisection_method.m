function [root,count] = bisection_method(f,a,b,e)
     count=0;
    while abs(b-a) >= e
        m = (a + b)/2;
        count = count + 1;
        if f(a) * f(m) <= 0
            b = m;
        else
            a = m;
        end
    end
    root = (a + b)/2;
end