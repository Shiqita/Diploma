function res=df1(xx,roo)
    syms f(x, ro)
    fx1=diff(f(x,ro),x,1);
    res=double(subs(fx1,{x,ro},{xx, roo}));
end