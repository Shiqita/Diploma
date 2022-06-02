%function mainfun()

%initial vars

xi=(2:17:1000);

%sobstvenniey zna4eniya
h=49; 
l=75; 

ro=sqrt(2); 
w=sqrt(1000); 

n=length(xi);

%eps=;

%% setup 
r=zeros(n);
teta=zeros(n);
omega=zeros(n);
L=zeros(n);
dL=zeros(n);
%-------------
arg=[0 ro w l h];
teta(1)=-v2(arg);
%-------------

r_=zeros(n);
teta_=zeros(n);
omega_=zeros(n);
L_=zeros(n);
dL_=zeros(n);
%-------------
teta_(n)=-pi/2;
%-------------

X=set_x(xi,ro); % poly4aem iz XI

%NormControl 'on' 'off' - mnogocomponentnoe vi4islene oshibki|refine n?
options=odeset('RelTol',eps,'Refine',10,'Stats','on');

%%
%TETA_OMEGA ->
tic

args=[ro h l w teta(1) 1];
sol=ode45(set_teta_omega(x,args),(X(1):X(2)),args,options);
temp=deval(sol,sol.x);
omega(1)=temp(2,end);
teta(1)=temp(1,end);
    
for i=2:n-1
    args=[ro h l w teta(i-1) 1];
    sol=ode45(set_teta_omega(x,args),(X(i):X(i+1)),args,options);
    temp=deval(sol,sol.x);
    omega(i)=temp(2,end);
    teta(i)=temp(1,end);
end
%r <-
r=R(omega,n,1);
disp('TETA+OMEGA+r')
toc
%L L' ->
tic
L=Lambda(x,ro,w,l,h,teta,r,1,n);
dL=dLambda(x,ro,w,l,h,teta,L,1,n);
toc

%%
%TETA_OMEGA <-
tic

args=[ro h l w teta_(n) 1];
sol_=ode45(set_teta_omega(x,args),(X(n):X(n-1)),args,options);
temp=deval(sol_,sol_.x);
omega_(1)=temp(2,end);
teta_(1)=temp(1,end);
    
for i=n-1:1
    args(5)=teta_(i+1);
    sol_=ode45(set_teta_omega(x,args),(X(i+1):X(i)),args,options);
    temp=deval(sol_,sol_.x);
    omega_(i)=temp(2,end);
    teta_(i)=temp(1,end);
end
%r ->
r_=R(omega_,1,n);
disp('TETA+OMEGA+r')
toc
%L L' <-
tic
L_=Lambda(x,ro,w,l,h,teta_,r_,n,1);
dL_=dLambda(x,ro,w,l,h,teta_,L_,n,1);
toc

% massiv must be filled with ~1
disp(check(r,r_,tet,tet_));



   % arg=[x ro w l h];  v funciyah
 


%% X|f|f'|f''                             real  functions
function res=set_x(xi,ro)
    res=zeros(length(xi));
    for i=1:length(xi)
        res(i)=sqrt(xi(i)-ro^2);
    end
end

function fx=f(x, ro)
    fx=((x^2+ro^2)*(x^2+ro^2-1))^(1/2);
end

function res=df1(xx,roo)
    syms f(x, ro)
    fx1=diff(f(x,ro),x,1);
    res=double(subs(fx1,{x,ro},{xx, roo}));
end

function res=df2(xx,roo)
    syms f(x, ro)
    fx2=diff(f(x,ro),x,2);
    res=double(subs(fx2,{x,ro},{xx, roo}));
end
%% q|U|psi podgotovka k TETA|OMEGA
function qx=q(x,ro,h,l,w)
    qx=h*ro^2-l*(ro^2)*(x^2+ro^2)+(w^2)*(x^2+ro^2)^2;
end

function ux=U(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    p1=q(x,h,l,ro,w)/(f(x,ro))^2;
    p2=-(df2(x,ro))/(2*f(x,ro));
    p3=(df1(x,ro)^2)/(4*f(x,ro)^2);
    ux=p1+p2+p3;
end

%peredaem w-omega o-teta(end)
function ps=psi(arg,teta)
    x=arg(1);
    w=arg(3);
    ps=w*x+v2(arg)+teta;
end
%% TETA|OMEGA
% TETA
function res=TETA(arg,args)
    tet=args(5);

    p1=((U(arg)/(v2(arg))^2-w-dv2(arg))*cos(psi(arg,tet))^2);
    p2=(v1(arg)^2-w-dv2(arg))*sin(psi(arg,tet))^2;
    p3=-((dv1(arg))/(2*v1(arg)^2))*sin(2*psi(arg,tet));
    res=p1+p2+p3;
end

function res=V(arg,tet)
    p1=(dv1(arg))/(v1(arg))*cos(2*psi(arg,tet));
    p2=((U(arg)/v1(arg)^2)-v1(arg)^2)*sin(2*psi(arg,tet));
    res=(p1+p2)/2;
end

function res=set_teta_omega(x,args)
    arg=[x args(1) args(4) args(3) args(2)];  
    om=args(end);
    args(end)=[];
    
    tet=TETA(arg,args);
    tmp=-V(arg,tet)*om;
    res=[tet; tmp];
end
%% R|f*Lambda|f'*Lambda' !! start-1 <-  !!!  start+1 -> !!
function res=R(Omega,start,endl)
    res(start)=1;
    if start>endl
        k=1; %start=n-1 endl=1
        for i=start-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    else
        k=-1; %start=2 endl=n
        for i=start-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    end
end

function res=Lambda(x,ro,w,l,h,teta,r,start,endl)
    res=zeros(length(x));
    arg=[0 ro w l h];
    for i=start:endl
        arg(1)=x(i);
        p1=(r(i)/v1(arg));
        p2=cos(w*x(i)+v2(arg)+teta(i));
        res(i)=(p1*p2)/(f(x(i),ro))^(1/2);
    end
end

function res=dLambda(x,ro,w,l,h,teta,L,start,endl)
    res=zeros(length(x));
    arg=[0 ro w l h];
    for i=start:endl
       arg(1)=x(i);
       p1=(L(i)*f(x(i),ro)^(1/2))/(df1(x(i),ro)^(1/2));
       p2=-(v2(arg)^2)*tan(w*x(i)+v2(arg)+teta(i));
       res(i)=p1*p2;
    end    
end
%% error check 
function res=check(r,r_,tet,tet_)
    res=zeros(length(r));
    for i=1:length(r)
        res(i)=r(i)*r_(i)*sin(tet(i)-tet_(i));
    end
end
%%                                             end of real functions
%v1|v2|dv1|dv2
function res=v1(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    P=1+P1/(x^2+ro^2)+P2/(x^2+ro^2)^2;
    res=sqrt(w*P);
end

function res=v2(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    p1=(P1/ro)*atan(x/ro);
    tmp=x/(x^2+ro^2)+(1/ro)*atan(x/ro);
    p2=(P2/(2*ro^2))*tmp;
    res=w*(p1+p2);
end

%dv1=(v1^2)'
function res=dv1(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);   
    
    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);

    p1=-(2*P1*x)/(x^2+ro^2)^2;
    p2=-(4*P2*x)/(x^2+ro^2)^3;

    res=w*(p1+p2);
end
function res=dv2(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);
    
    P1=R_1(w,ro,l);
    P2=R_2_(w,ro,l,h);
    
    tmp=(1/(ro^2*((x^2)/(ro^2)+1))+(ro^2-x^2)/(x^2+ro^2)^2);
    p2=(P2/(2*ro^2))*tmp;
    p3=P1/(ro^2*((x^2)/(ro^2)+1));
    
    res=w*(p2+p3);
end
%% R1|R2|R2_
%                eto L!
function res=R_1(w,ro,l)
    res=(w^2-l*ro^2)/(2*w^2);
end
function res=R_2(w,ro,l,h)
%         eto L!
    p=R_1(w,ro,l);
    P1=(h*ro^2)/(2*w^2);
    P2=p*(1-ro^2-p/2);
    P3=-(2*ro^2-1)/(4*w^2);
    res=P1+P2+P3;
end
function res=R_2_(w,ro,l,h)
    res=R_2(w,ro,l,h)+R_1(w,ro,l)*ro^2;
end
%%                                              