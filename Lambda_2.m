function Lambda_2(sqr_ro,sqr_w,l,h,last_tet,step,end_xi)
    %sobstvenniey zna4eniya
    h=575.34;
    l=1228.06;

    ro=sqrt(2);
    w=sqrt(1000);

    xi=(ro^2:0.01:18);

    X=set_x(xi,ro);
    n=length(X);
    %% setup 
    temp_R=set_R(w,ro,l,h); %Set R1 R2 ^R2
    P1=temp_R(1);
    P2=temp_R(2);
    P3=temp_R(3);
    %-------------
    teta=zeros(n,1);
    omega=zeros(n,1);
    omega(n)=1;
    teta(n)=last_tet-pi/2;                                
    %-------------

    options=odeset('AbsTol',1e-06,'RelTol',1e-03,'Refine','8');
    %%
    %TETA_OMEGA <-
    for i=n-1:-1:1
        y0=[teta(i+1) 1];
        tspan=[X(i+1) X(i)];
        disp(tspan);

        sol=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
        temp=deval(sol,sol.x);
        omega(i)=temp(2,end);
        teta(i)=temp(1,end);
    end

    %r ->
    r=R(omega,1,n);
    %L L' <-
    L=Lambda(X,ro,w,l,h,P1,P2,P3,teta,r,n,1);
    dL=dLambda(X,ro,w,l,h,P1,P2,P3,teta,L,n,1);

    plot(xi,L);
    hold;
end
%plot(X,dL)
%%                                              main  functions
%% f|f^2|f'|d(sqrt(f))|f''                          
function res=f(x, ro)
    res=sqrt((x^2+ro^2)*(x^2+ro^2-1));
end
function res=f_1(x, ro)
    res=(x^2+ro^2)*(x^2+ro^2-1);
end
function res=df1(x, ro)
    res=(2*x*(ro^2 + x^2) + 2*x*(ro^2 + x^2 - 1))/(2*((ro^2 + x^2)*(ro^2 + x^2 - 1))^(1/2));
end
function res=df12(x, ro)
    res=(2*x*(ro^2 + x^2) + 2*x*(ro^2 + x^2 - 1))/(4*((ro^2 + x^2)*(ro^2 + x^2 - 1))^(3/4));
end
function res=df2(x, ro)
    a=(4*ro^2 + 12*x^2-2)/(2*((ro^2 + x^2)*(ro^2+x^2-1))^(1/2));
    b=-(2*x*(ro^2+x^2)+2*x*(ro^2+x^2-1))^2/(4*((ro^2+x^2)*(ro^2+x^2-1))^(3/2));
    res=a+b;
end
%% TETA|OMEGA
function res=set_teta_omega(x,y,ro,w,l,h,P1,P2,P3)   
    arg=[x ro w l h P1 P2 P3];
%v1|v2----------------------------------
    V1=(v1(arg));
    V2=(v2(arg));
%dv1|dv2--------------------------------
    DV1=dv1(arg);    
    DV2=dv2(arg);
%psi------------------------------------
    psi=w*x+V2+y(1);
%U--------------------------------------
    K0=f(x,ro);
    K1=f_1(x,ro);
    K2=df1(x,ro);
    K3=df2(x,ro);
    
    q=(h*ro^2)-(l*(ro^2))*(x^2+ro^2)+(w^2)*(x^2+ro^2)^2;
    
    U=(q/K1)-(K3/(2*K0))+((K2^2)/(4*K1));
%TETA-----------------------------------    
    t11=roundn(((U/V1)-w-DV2),-9);
    t12=roundn((cos(psi))^2,-9);
    t1=t11*t12;
    
    t21=roundn((V1-w-DV2),-9);
    t22=roundn(sin(psi)^2,-9);
    t2=t21*t22;
    
    t31=-roundn((DV1/(2*V1)),-9);
    t32=roundn(sin(2*psi),-9);
    t3=t31*t32;
%---------------------------------------    
    F1=t1+t2+t3;
%V--------------------------------------        
    Om11=roundn((DV1/V1),-9);
    Om12=roundn(cos(2*psi),-9);
    Om1=Om11*Om12;
    
    Om21=roundn(((U/V1)-V1),-9);
    Om22=roundn(sin(2*psi),-9);
    Om2=Om21*Om22;
%---------------------------------------    
    F2=-((Om1+Om2)/2)*y(2);
%---------------------------------------
    res=[F1; F2]; 
end
%% R|f*Lambda|f'*Lambda' !! start-1 <-  !!!  start+1 -> !!
function res=R(Omega,start,endl)
    res(start)=1;
    if start>endl
        k=1; %start=n-1 endl=1
        for i=start-k:-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    else
        k=-1; %start=2 endl=n
        for i=start-k:-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    end
end

function res=Lambda(x,ro,w,l,h,P1,P2,P3,tet,rr,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h P1 P2 P3];
    if start>endl
        k=-1;
    else
        k=1;
    end
    for i=start:k:endl
        arg(1)=x(i);
        psi=w*x(i)+v2(arg)+tet(i);
        
        temp=(rr(i)/sqrt(v1(arg)))*cos(psi);        %edinstvennoe mesto gde v1 ne v ^2
        res(i)=temp/sqrt(f(x(i),ro));
    end
end

function res=dLambda(x,ro,w,l,h,P1,P2,P3,tet,LL,start,endl)
    res=zeros(length(x),1);
    arg=[0 ro w l h P1 P2 P3];
    if start>endl
        k=-1;
    else
        k=1;
    end
    for i=start:k:endl
       arg(1)=x(i);
       psi=w*x(i)+v2(arg)+tet(i);
       df=df12(x(i),ro);
       
       p1=(LL(i)*sqrt(f(x(i),ro)))/(df);
       p2=-(v1(arg))*tan(psi);
       res(i)=p1*p2;
    end    
end
%%                                       end of real functions
%v1|v2
function res=v1(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);
    res=w*(1+P1/(x^2+ro^2)+P3/((x^2+ro^2)^2));
end
function res=v2(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);
    res=w*((P1/ro)*atan(x/ro)+(P3/(2*ro^2))*(x/(x^2+ro^2)+(atan(x/ro)/ro)));
end

%%  d(v1^2)|d(v2)
%dv1=(v1^2)'
function res=dv1(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);   
    res=-w*((2*P1*x)/(ro^2+x^2)^2+(4*P3*x)/(ro^2+x^2)^3);
end
function res=dv2(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);
    res=w*(P1/(ro^2*(x^2/ro^2+1))+(P3*(1/(ro^2+x^2)-(2*x^2)/(ro^2+x^2)^2+1/(ro^2*(x^2/ro^2+1))))/(2*ro^2));
end
%% R1|R2|R2_
%                eto  L!
function res=R_1(w,ro,l)
    res=(w^2-(l*ro^2))/(2*w^2);
end
function res=R_2(w,ro,h,p)
    P1=(h*ro^2)/(2*w^2);
    P2=p*(1-ro^2-p/2);
    P3=-((2*ro^2-1)/(4*w^2));
    res=P1+P2+P3;
end
function res=R_2_(ro,p1,p2)
    res=p2+p1*ro^2;
end   
function res=set_R(w,ro,l,h)
    P1=R_1(w,ro,l);
    P2=R_2(w,ro,h,P1);
    P3=R_2_(ro,P1,P2);
    
    res=[P1 P2 P3];
end
function res=set_x(xi,ro)
    res=zeros(length(xi),1);
    for i=1:length(xi)
        res(i)=sqrt(xi(i)-ro^2);
    end
end
%%                                             end of main functions