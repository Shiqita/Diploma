function res=dLambda(x,ro,w,l,h,teta,L,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h];
    for i=start:endl
        arg(1)=x(i);
        p1=(L(i)*f(x(i),ro)^(1/2))/(df1(x(i),ro)^(1/2));
        p2=-(v2(arg)^2)*tan(w*x(i)+v2(arg)+teta(i));
        res(i)=p1*p2;
    end    
end

