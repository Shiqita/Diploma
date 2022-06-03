function res=Lambda(x,ro,w,l,h,teta,r,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h];
    for i=start:endl
        arg(1)=x(i);
        p1=(r(i)/v1(arg));
        p2=cos(w*x(i)+v2(arg)+teta(i));
        res(i)=(p1*p2)/(f(x(i),ro))^(1/2);
    end
end

