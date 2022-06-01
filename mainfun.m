%function mainfun()

%initial vars

xi=(2:17:1000);

%h=; 
%l=; 

ro=sqrt(2); 
%w=; 

%n=;    obshee kol-vo n v x0-x1

%eps=;

%% setup

teta=zeros(length(xi));
teta(1)=-v2();

X=set_x(xi,ro); % poly4aem iz XI

 
%NormControl 'on' 'off' - mnogocomponentnoe vi4islene oshibki
options=odeset('RelTol',eps,'Refine',n,'Stats','on');

%%
%New TETA and OMEGA after 30.05 Omega at start always 1
for i=1:n
    args=[ro h l w teta(i) 1];
    sol=ode45(@() set_teta_omega(x,args),(X(i):X(i+1)),args,options);
    temp=deval(sol,sol.x);
    omega(i)=temp.2;
    teta(i)=temp.1;
end
%% Omega and TETA before 30.05

tic
%TETA
[t0,teta]=ode45(@(t0,teta) TETA(x, args),(t0:tn),[ro h l w teta(end)]);
plot(t0,teta,'-b');

toc

tic

%OMEGA
Omega0=1;
Omega=zeros(length(t0),length(t0));
for i=1:length(t0)
[t1,Omega(i)]=ode45(@(t1,Omega(i)) OMEGA(x,args),(),[ro h l w teta() Omega0]);

% ???
end
toc
%end


%      test functions

%      test functions
%%    real  functions
%%
%fffff
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

%fffff
%%
function qx=q(x,ro,h,l,w)
qx=h*ro^2-l*(ro^2)*(x^2+ro^2)+(w^2)*(x^2+ro^2)^2;
end

function ux=U(x,ro,h,l,w)
ux=(q(x,h,l,ro,w)/(f(x,ro))^2-(df2(x,ro))/(2*f(x,ro))+(df1(x,ro)^2)/(4*f(x,ro)^2));
end

%peredaem w-omega o-teta(end)
function ps=psi(x,w,teta)
ps=w*x+v2(x)+teta;
end
%%
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

function res=V(x,teta)
%p1=([v1(x)^2)]'/(v1(x))*cos(2*psi(x,w,teta));
% [v1(x)^2]'   !!!!!!!!
p2=((U()/v1(x)^2)-v1(x)^2)*sin(2*psi(x,w,teta));
res=(p1+p2)/2;
end

% OMEGA            TETA moget bit' zamenena na var
function res=OMEGA(x,args)
om=args(end);
args(end)=[];
res=-V(x,TETA(x,args))*om;
end

function res=set_teta_omega()

end
%%
function res=R(x,args)
r=args(end);
args(end)=[];
res=V(x,TETA(x,args))*r;
end
%% end of real functions  