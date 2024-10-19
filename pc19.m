function [u,t,counter]=pc19(f,t0,y0,tf,tol,maxord)
%   Solve non-stiff differential equations, variable order-variable step method.
%   
%   [UOUT,TOUT,COUNTER]=pc19(ODEFUN,T0,Y0,TFINAL,TOL) integrates the system of
%   differential equations y'=f(t,y) from time T0 to TFINAL with initial
%   conditions Y0 and an error tolerance TOL. ODEFUN is a function handle. 
%   For a scalar T and vector Y, ODEFUN(T,Y) must return a column vector 
%   corresponding to f(t,y). Each row in UOUT corresponds to a time 
%   returned in the column vector TOUT. Each column in UOUT corresponds 
%   to a solution of one differential equation in the system f(t,y).
%   COUNTER specifies the number of function evaluations.
%
%   [UOUT,TOUT]=pc19(ODEFUN,T0,Y0,TFINAL,TOL,MAXORD) solves as above but
%   limits the maximum order of the method to MAXORD. By default,
%   MAXORD=9.
%   
%   This function utilizes an adaptive-step scheme based on Adams-Bashforth
%   and Adams-Moulten methods for a predictor-corrector strategy. The
%   initial step utilizes Heun's Method. The proceeding kth step utilizes a
%   kth order method, beginning at 2nd order and maximizing at 9th order by
%   default. In the case MAXORD is inputted, the order maximizes at MAXORD.
%   In either case, the order remains at MAXORD (9) until TFINAL is
%   reached.
%

addpath(genpath(pwd));

if nargin<6
    maxord=9;
end

% Initialize function calls
getParams = {@(h,H) getParams2(h,H), @(h,H) getParams3(h,H), ...
             @(h,H) getParams4(h,H), @(h,H) getParams5(h,H), ...
             @(h,H) getParams6(h,H), @(h,H) getParams7(h,H), ...
             @(h,H) getParams8(h,H), @(h,H) getParams9(h,H)};

% INITALIZATION
t(1)=t0;
u=y0;
H=(tf-t0)/100; % Choose an initial step size
tol1=tol/(tf-t0); % Bound on the LTE (error committed on one step)
etarg=tol1/2; % Target error per step is within half of tol1.
counter=0;
[u,t,km,counter]=OneStep(f,u,t,tol1,H,counter); % Take first step with Heun's Method
H=t(2)-t(1); % Get latest H

% MAIN ITERATION:
hmin=min(tol/100,1e-6);
j=2;
goodstep=true;
while t(j)<tf && H>hmin
    % Get the parameters of the algorithm, based on the current and
    % previous step sizes (H, h(1), ... , h(k-1) respectively):
    ord=min(j,maxord); % Order increases with each step up until maxord
    h=getSteps(ord,t,j); % Get step sizes
    [pbeta,cbeta,err]=getParams{ord - 1}(h, H); % Get parameters
    if goodstep
        km(:,ord)=f(t(j),u(:,j)); % Compute most recent slope
        counter=counter+1;
    end
    u1=u(:,j)+km*pbeta; % Compute Predictor
    slp=[km(:,2:ord),f(t(j)+H,u1)]; % Shift slopes to insert implicit slope.
    counter=counter+1;
    u2=u(:,j)+slp*cbeta; % Compute Corrector
    plte=err(2)*(u2-u1)/(H*(err(1)-err(2))); % Estimate PLTE
    % Decide how to proceed (step-size cannot change by a factor > 10):
    if norm(plte)>tol1
        % Reject this step, reduce H, and try again:
        goodstep=false;
        ht=estimateNextH(err,ord,u1,u2,etarg,h);
        H=max(ht,0.1*H);
        
    else
        % Accept step
        goodstep=true;
        t(j+1)=t(j)+H; % Take step forward
        u(:,j+1)=u2+H*plte; % Local Extrapolation
        j=j+1;
        h=getSteps(ord,t,j);
        ht=min(estimateNextH(err,ord,u1,u2,etarg,h),tf-t(j));
        if ht>10*H
            H=10*H;
        elseif ht<0.1*H
            H=0.1*H;
        else
            H=ht;
        end
        if t(j)<tf && ord==maxord % Shift arrays to accomdate next step
            km=km(:,2:ord);
        end
    end
    %fprintf('t=%.8e, h=%.4e\n',t(j),H);
    if H<=hmin && tf-t(j)>H
        fprintf('Minimum step size reached; integration terminated.\n')
    end 
end
t=t';
u=u';
end