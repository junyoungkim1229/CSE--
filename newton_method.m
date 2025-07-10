function [root,count]=newton_method(f,df,x0,e,iMax)

for i=1:iMax
    x=x0-f(x0)/df(x0);%Newton method
    if abs((x-x0))<e  %error
        root=x; count=i;
        break
    end
    x0=x;

    if i==iMax
        disp("error")
    end
   
end
