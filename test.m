GLOBAL args

args=[1 2 3 4 5 6 7];
disp(args)

X=[1 1 1 1 1 1 1];

sol=ode45(@(x,y) set_teta_omega(x,args),1:10,args);


function res=set_teta_omega(x,args)
    tmp1=0;
    tmp2=0;
    for i=1:length(x)-1
        tmp1=+x(i);
        tmp2=+args(i+1);
    end
    res=[tmp1;tmp2];
end