%function AdaptiveRunge(h,n,y0)
    clc;clear;close all;
    % initial conditions
    t=0;y=0;
    linet=(0);liney=(0);
    t0=t;tn=t0+4;y0=y;y1=1.2;
    % step size=%f\n
    h=0.01;
    % tolerance
    dt_min=0.001;
    y_tol=0.001;
    dy_min=0.008;
    dy_max=0.01;
    % start response plot
    figure()
    plot(t,y);hold;
    tic
    while t<tn
        % current step
        % calculate partial steps
        k1=mfun(t,y);
        k2=mfun(t+h/2,y+h*k1/2);
        k3=mfun(t+h/2,y+h*k2/2);
        k4=mfun(t+h,y+h*k3);
        % combine partial steps
        sy=y+h/6*(k1+2*k2+2*k3+k4);

        % half step
        % calculate partial steps
        k2=mfun(t+h/4,y+h*k1/4);
        k3=mfun(t+h/4,y+h*k2/4);
        k4=mfun(t+h/2,y+h*k3/2);
        % combine partial steps
        hsy=y+h/12*(k1+2*k2+2*k3+k4);
        
        % double step
        % calculate partial steps
        k2=mfun(t+h,y+h*k1);
        k3=mfun(t+h,y+h*k2);
        k4=mfun(t+h*2,y+h*k3*2);
        % combine partial steps
        dsy=y+h/3*(k1+2*k2+2*k3+k4);

        %     % compare and use
        if abs(sy)<y_tol
            if h~=dt_min
                fprintf('New step size=%f\n',dt_min);
                h=dt_min;
            end
            ny=sy;
        else
            if abs(sy)>y_tol && (abs(sy-hsy)/abs(sy))>dy_max
                h=h/2;
                fprintf('New step size=%f\n',h);
                ny=hsy;
            elseif abs(sy)>y_tol && (abs(sy-dsy)/abs(sy))<dy_min
                h=h*2;
                fprintf('New step size=%f\n',h);
                ny=dsy;
            else
                ny=sy;
            end
        end
        y=ny;
        t=t+h;
        linet=[linet t+h];
        liney=[liney ny];
    end
  % plot(linet,liney)
    toc
    tic
    solution=ode45(@mfun,(t0:h:tn),[y0 y1]);
    disp(solution.x)
    plot(solution.x,solution.y)
    toc
    function mf=mfun(t,y)
    mf=exp(y)-y*exp(t);
    %mf=sqrt(tan(t)+y);
    %mf=sin(y)+cos(sqrt(t)*y);
    end
%end
