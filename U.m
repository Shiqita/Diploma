function ux=U(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    p1=q(x,h,l,ro,w)/(f(x,ro))^2;
    p2=-(df2(x,ro))/(2*f(x,ro));
    p3=(df1(x,ro)^2)/(4*f(x,ro)^2);
    
    ux=p1+p2+p3;
end