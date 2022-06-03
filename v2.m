function res=v2(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    
    p1=(P1/ro)*atan(x/ro);
    tmp=x/(x^2+ro^2)+(1/ro)*atan(x/ro);
    p2=(P2/(2*ro^2))*tmp;
    
    res=w*(p1+p2);
end