clear all; close all; clc;
y0=[1.2;0;0;-1.04935751];
t0=0;
tf=2*pi;
tol=1e-6;
tol1=2e-5*tol;
tspan=[t0 tf];
opt=odeset('AbsTol',tol1,'RelTol',tol1);

r13=@(t,y)((y(1)+1/82.45)^2+y(2)^2).^1.5;
r23=@(t,y)((y(1)-(1-1/82.45))^2+y(2)^2).^1.5;
f=@(t,y)[y(3);
    y(4);
    2*y(4)+y(1)-(1-1/82.45)*((y(1)+1/82.45)/r13(t,y))-1/82.45*((y(1)-(1-1/82.45))/r23(t,y));
    -2*y(3)+y(2)-(1-1/82.45)*(y(2)/r13(t,y))-1/82.45*(y(2)/r23(t,y))];

[u,t,counter]=pc113(f,t0,y0,tf,tol,9);
SOL=ode45(f,tspan,y0,opt);
SOL1=ode113(f,tspan,y0,opt);
t1=SOL.x';
u1=SOL.y';
t2=SOL.x';
u2=SOL.y';
odecounts=SOL.stats.nfevals;
odecounts1=SOL1.stats.nfevals;
fprintf('pc113 function evals: %f\n',counter)
fprintf('ode45 function evals: %f\n',odecounts)
fprintf('ode113 function evals: %f\n',odecounts1)

figure(1);
plot(u(:,1),u(:,2),Color='m',Marker='diamond'); hold on; grid on;
plot(u1(:,1),u1(:,2),Color='c',Marker='+');
plot(u2(:,1),u2(:,2),Color='y',Marker='o');
title('Orbital Problem Solution');
ylabel('y(t)');
xlabel('x(t)');

