function res=V(arg,tet)
    p1=(dv1(arg))/(v1(arg))*cos(2*PSI(arg,tet));
    p2=((U(arg)/v1(arg)^2)-v1(arg)^2)*sin(2*PSI(arg,tet));
    res=(p1+p2)/2;
end

