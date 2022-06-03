function res=TETA(arg)
    w=arg(3);tet=arg(end);arg(end)=[];

    p1=((U(arg)/(v2(arg))^2-w-dv2(arg))*cos(psi(arg,tet))^2);
    p2=(v1(arg)^2-w-dv2(arg))*sin(psi(arg,tet))^2;
    p3=-((dv1(arg))/(2*v1(arg)^2))*sin(2*psi(arg,tet));
    res=p1+p2+p3;
end