function res=dv2(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);
    
    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    
    tmp=(1/(ro^2*((x^2)/(ro^2)+1))+(ro^2-x^2)/(x^2+ro^2)^2);
    p2=(P2/(2*ro^2))*tmp;
    p3=P1/(ro^2*((x^2)/(ro^2)+1));
    
    res=w*(p2+p3);
end