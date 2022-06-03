%dv1=(v1^2)'
function res=dv1(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);   
    
    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);

    p1=-(2*P1*x)/(x^2+ro^2)^2;
    p2=-(4*P2*x)/(x^2+ro^2)^3;

    res=w*(p1+p2);
end