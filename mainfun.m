%function mainfun()
%initial vars
%global h l w ro
xi=(2:250:35000);

%sobstvenniey zna4eniya
h=5; 
l=12; 

ro=sqrt(2); 
w=sqrt(10); 

n=length(xi);

%eps=;

%% setup 
%r=zeros(n,1);
teta=zeros(n,1);
omega=zeros(n,1);
%L=zeros(n,1);
%dL=zeros(n,1);
%-------------
arg=[0 ro w l h];
teta(1)=pi/2-v2(arg);                                          
%-------------

%r_=zeros(n,1);
teta_=zeros(n,1);
omega_=zeros(n,1);
%L_=zeros(n,1);
%dL_=zeros(n,1);
%-------------
teta_(n)=-pi/2;
%-------------

X=set_x(xi,ro); % poly4aem iz XI

%NormControl 'on' 'off' - mnogocomponentnoe vi4islene oshibki|refine n? 'RelTol',eps,
options=odeset('Refine',4,'Stats','on');

%%
%TETA_OMEGA ->
tic

%args=[ro h l w teta(1) 1];
toc
tic
sol=ode45(@(x,y) set_teta_omega(x,y,ro,w,h,l),(X(1):X(2)),[teta(1) 1],options);
toc
temp=deval(sol,sol.x);
omega(1)=temp(2,end);
teta(1)=temp(1,end);
    

for i=2:n-1
    tic
    disp(i)
    %args=[ro h l w teta(i-1) 1];
    sol=ode45(@(x,y) set_teta_omega(x,y,ro,w,h,l),(X(i):X(i+1)),[teta(i-1) 1],options);
    temp=deval(sol,sol.x);
    omega(i)=temp(2,end);
    teta(i)=temp(1,end);
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

sol_=ode45(@(x,y) set_teta_omega(x,y,ro,w,h,l),(X(n):X(n-1)),[teta_(n) 1],options);
temp=deval(sol_,sol_.x);
omega_(n)=temp(2,end);
teta_(n)=temp(1,end);

for i=n-1:1
    sol_=ode45(@(x,y) set_teta_omega(x,y,ro,w,h,l),(X(i+1):X(i)),[teta_(i+1) 1],options);
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
L_=Lambda(X,ro,w,l,h,teta_,r_,n,1);
dL_=dLambda(X,ro,w,l,h,teta_,L_,n,1);
toc

% massiv must be filled with ~1
disp(check(r,r_,tet,tet_));

%% X|f|f'|f''                             real  functions

function res=set_teta_omega(x,y,ro,w,h,l)
    arg=[x ro w l h];
    om=y(2);
    
    tet=TETA([arg y(1)]);
    tmp=-V(arg,tet)*om;
    res=[tet; tmp];
end
%% R|f*Lambda|f'*Lambda'

%% error check 
function res=check(r,r_,tet,tet_)
    res=zeros(length(r),1);
    for i=1:length(r)
        res(i)=r(i)*r_(i)*sin(tet(i)-tet_(i));
    end
end
%%                                             end of real functions                             