function res=v1(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    P=1+P1/(x^2+ro^2)+P2/(x^2+ro^2)^2;
    res=sqrt(w*P);
end