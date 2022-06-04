%function mainfun()
%initial vars

xi=(0:0.5:10);

%sobstvenniey zna4eniya
h=5; 
l=12; 

ro=sqrt(2); 
w=sqrt(6); 

n=length(xi);

%eps=;

%% setup 
temp_R=set_R(w,ro,l,h); %Set R1 R2 ^R2
arg=[0 ro w l h temp_R];
P1=temp_R(1);
P2=temp_R(2);
P3=temp_R(3);
%-------------
teta=zeros(n,1);
omega=zeros(n,1);
omega(1)=1;
teta(1)=pi/2-v2(arg);                                          
%-------------
teta_=zeros(n,1);
omega_=zeros(n,1);
omega_(n)=1;
teta_(n)=-pi/2;
%-------------

X=set_x(xi,ro); % poly4aem iz XI

%NormControl 'on' 'off' - mnogocomponentnoe vi4islene oshibki|refine n?
%'RelTol',eps, 'Refine',4,
options=odeset('Stats','on');

%%
%TETA_OMEGA ->
for i=1:n-1
    tic
    disp(i)
    y0=[teta(i) 1];
    tspan=[X(i) X(i+1)];
    disp(tspan);
    
    sol=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
    temp=deval(sol,sol.x);
    omega(i+1)=temp(2,end);
    teta(i+1)=temp(1,end);
    %[x,y]=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
    %omega(i+1)=y(2,end);
    %teta(i+1)=y(1,end);
    disp(omega(i+1));
    disp(teta(i+1));
    toc
end

%r <-
r=R(omega,n,1);
disp('TETA+OMEGA+r')
toc

%L L' ->
tic
L=Lambda(X,ro,w,l,h,teta,r,1,n);
dL=dLambda(X,ro,w,l,h,teta,L,1,n);
toc

%%
%TETA_OMEGA <-
tic
for i=n-1:1
    tic
    disp(i)
    y0=[teta(i+1) 1];
    tspan=[X(i+1) X(i)];
    disp(tspan);
    
    sol_=ode45(@(x,y) set_teta_omega(x,y,ro,w,l,h,P1,P2,P3),tspan,y0,options);
    temp=deval(sol_,sol_.x);
    omega_(i)=temp(2,end);
    teta_(i)=temp(1,end);
    disp(i)
    disp(omega);
    disp(teta);
    toc
end

%r ->
r_=R(omega_,1,n);
disp('TETA+OMEGA+r')
toc

%L L' <-
tic
L_=Lambda(X,ro,w,l,h,teta_,r_,n,1);
dL_=dLambda(X,ro,w,l,h,teta_,L_,n,1);
toc

% massiv must be filled with ~1
disp(check(r,r_,tet,tet_));

%%                                              main  functions
%% X|f|f'|f''                             
function res=set_x(xi,ro)
    res=zeros(length(xi),1);
    for i=1:length(xi)
        res(i)=abs(sqrt(xi(i)-ro^2));
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
function res=q(x,ro,h,l,w)
    res=h*ro^2-l*(ro^2)*(x^2+ro^2)+(w^2)*(x^2+ro^2)^2;
end

function res=U(arg)
    x=arg(1);ro=arg(2);w=arg(3);l=arg(4);h=arg(5);

    p1=q(x,h,l,ro,w)/(f(x,ro))^2;
    p2=-(df2(x,ro))/(2*f(x,ro));
    p3=(df1(x,ro)^2)/(4*f(x,ro)^2);
    res=p1+p2+p3;
end

function res=psi(arg,teta)
    x=arg(1);w=arg(3);
    res=w*x+v2(arg)+teta;
end
%% TETA|OMEGA
% TETA
function res=TETA(arg)
    w=arg(3);tet=arg(end);arg(end)=[];

    p1=((U(arg)/(v2(arg))^2-w-dv2(arg))*cos(psi(arg,tet))^2);
    p2=(v1(arg)^2-w-dv2(arg))*sin(psi(arg,tet))^2;
    p3=-((dv1(arg))/(2*v1(arg)^2))*sin(2*psi(arg,tet));
    res=p1+p2+p3;
end

function res=V(arg)
    tet=arg(end);arg(end)=[];
    
    p1=(dv1(arg))/(v1(arg))*cos(2*psi(arg,tet));
    p2=((U(arg)/v1(arg)^2)-v1(arg)^2)*sin(2*psi(arg,tet));
    res=(p1+p2)/2;
end

function res=set_teta_omega(x,y,ro,w,l,h,P1,P2,P3)
    %tic;
    omega0=y(2); 
    arg=[x ro w l h P1 P2 P3];  
    
    teta=TETA([arg y(1)]);
    omega=-V([arg teta])*omega0;
    res=[teta; omega]; 
    %show=[x teta omega];
    %disp('x(i)  teta   omega');
    %disp(show);
    %toc;
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
    res=zeros(length(x),1);
    arg=[0 ro w l h];
    for i=start:endl
        arg(1)=x(i);
        p1=(r(i)/v1(arg));
        p2=cos(w*x(i)+v2(arg)+teta(i));
        res(i)=(p1*p2)/(f(x(i),ro))^(1/2);
    end
end

function res=dLambda(x,ro,w,l,h,teta,L,start,endl)
    res=zeros(length(x),1);
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
    res=zeros(length(r),1);
    for i=1:length(r)
        res(i)=r(i)*r_(i)*sin(tet(i)-tet_(i));
    end
end
%%                                             end of real functions
%v1|v2|dv1|dv2
function res=v1(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);
    
    P=1+P1/(x^2+ro^2)+P3/(x^2+ro^2)^2;
    res=sqrt(w*P);
end

function res=v2(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);

    p1=(P1/ro)*atan(x/ro);
    tmp=x/(x^2+ro^2)+(1/ro)*atan(x/ro);
    p2=(P3/(2*ro^2))*tmp;
    res=w*(p1+p2);
end

%dv1=(v1^2)'
function res=dv1(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);   

    p1=-(2*P1*x)/(x^2+ro^2)^2;
    p2=-(4*P3*x)/(x^2+ro^2)^3;

    res=w*(p1+p2);
end
function res=dv2(arg)
    x=arg(1);ro=arg(2);w=arg(3);P1=arg(6);P3=arg(8);
    
    tmp=(1/(ro^2*((x^2)/(ro^2)+1))+(ro^2-x^2)/(x^2+ro^2)^2);
    p2=(P3/(2*ro^2))*tmp;
    p3=P1/(ro^2*((x^2)/(ro^2)+1));
    
    res=w*(p2+p3);
end
%% R1|R2|R2_
%                eto L!
function res=R_1(w,ro,l)
    res=(w^2-l*ro^2)/(2*w^2);
end
function res=R_2(w,ro,h,p)
    P1=(h*ro^2)/(2*w^2);
    P2=p*(1-ro^2-p/2);
    P3=-(2*ro^2-1)/(4*w^2);
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
%%                                             end of main functions                             