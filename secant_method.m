function [root,count]=secant_method(f,x1,x2,e,iMax)

for i=1:iMax
    Xi=x2-f(x2)*(x1-x2)/(f(x1)-f(x2));
    if abs((Xi-x2))<e
        root=Xi; count=i;
        break
    end
    x1=x2;
    x2=Xi;
end
