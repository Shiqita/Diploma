%function mainfun()

%initial vars

xi=(2:17:1000);

%sobstvenniey zna4eniya
%h=; 
%l=; 

ro=sqrt(2); 
%w=; 

n=lenght(xi);    %obshee kol-vo n v x0-x1

%eps=;

%% setup
r=(length(xi));
teta=zeros(length(xi));
omega=zeros(length(xi));
L=zeros(length(xi));
dL=zeros(length(xi));

teta(1)=-v2(0);

X=set_x(xi,ro); % poly4aem iz XI

%NormControl 'on' 'off' - mnogocomponentnoe vi4islene oshibki|refine n?
options=odeset('RelTol',eps,'Refine',10,'Stats','on');

%%
%New TETA and OMEGA after 30.05 Omega at start always 1
tic

args=[ro h l w teta(1) 1];
sol=ode45(set_teta_omega(x,args),(X(1):X(2)),args,options);
temp=deval(sol,sol.x);
omega(1)=temp(2,end);
teta(1)=temp(1,end);
    
for i=2:n
    args=[ro h l w teta(i-1) 1];
    sol=ode45(@() set_teta_omega(x,args),(X(i):X(i+1)),args,options);
    temp=deval(sol,sol.x);
    omega(i)=temp(2,end);
    teta(i)=temp(1,end);
end
%r
r=R(omega);

disp('TETA+OMEGA+r')
toc

tic

L=Lambda();

dL=dLambda();

toc
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

function ux=U(x,ro,h,l,w)
    p1=q(x,h,l,ro,w)/(f(x,ro))^2;
    p2=-(df2(x,ro))/(2*f(x,ro));
    p3=(df1(x,ro)^2)/(4*f(x,ro)^2);
    ux=p1+p2+p3;
end

%peredaem w-omega o-teta(end)
function ps=psi(x,w,teta)
    ps=w*x+v2(x)+teta;
end
%% TETA|OMEGA
% TETA
function res=TETA(x,args)
    ro=args(1);
    h=args(2);
    l=args(3);
    w=args(4);
    tet=args(5);

    p1=((U(x,ro,h,l,w)/(v2(x))^2-w-dv2(x))*cos(psi(x,w,tet))^2);
    p2=(v1(x)^2-w-dv2(x))*sin(psi(x,w,tet))^2;
    % ??? p3=-((v1(x)^2)/(2*v1(x)^2))sin(2*psi(x,w,tet)) sho delat's
    % [v1(x)^2]'   !!!!!!!
    res=p1+p2+p3;
end

function res=V(x,tet)
    %p1=([v1(x)^2)]'/(v1(x))*cos(2*psi(x,w,tet));
    % [v1(x)^2]'   !!!!!!!!
    p2=((U()/v1(x)^2)-v1(x)^2)*sin(2*psi(x,w,tet));
    res=(p1+p2)/2;
end

function res=set_teta_omega(x,args)
    om=args(end);
    args(end)=[];
    tet=TETA(x,args);
    tmp=-V(x,tet)*om;
    res=[tet tmp];
end
%% R|f*Lambda|f*Lambda'
function res=R(Omega)
    res=zeros(length(Omega));
    res(end)=1;
    for i=(length(Omega)-1):1
        res(i)=Omega(i)*r;
    end
end

function res=Lambda(r,x,w,ro,teta)
    res=zeros(length(x));
    for i=1:length(x)
        p1=(r(i)/v1(x(i)));
        p2=cos(w*x(i)+v2(i)+teta(i));
        res(i)=(p1*p2)/(f(x(i),ro))^(1/2);
    end
end

function res=dLambda(x,w,ro,teta,L)
    res=zeros(length(x));
    for i=1:length(x)
       p1=(L(i)*f(x(i),ro))/df1(x(i),ro);
       p2=-v2(x(i))*tan(w*x(i)+v2(x(i)+teta(i)));
       res(i)=p1*p2;
    end    
end
%%                                                   end of real functions  