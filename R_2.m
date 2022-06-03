function res=R_2(w,ro,l,h)
%          eto L!
    p=R_1(w,ro,l);
    P1=(h*ro^2)/(2*w^2);
    P2=p*(1-ro^2-p/2);
    P3=-(2*ro^2-1)/(4*w^2);
    res=P1+P2+P3;
end