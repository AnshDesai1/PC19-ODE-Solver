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

[u,t,counter]=pc19(f,t0,y0,tf,tol,9);
SOL=ode45(f,tspan,y0,opt);
SOL1=ode113(f,tspan,y0,opt);
t1=SOL.x';
u1=SOL.y';
t2=SOL.x';
u2=SOL.y';
odecounts=SOL.stats.nfevals;
odecounts1=SOL1.stats.nfevals;
fprintf('pc19 function evals: %f\n',counter)
fprintf('ode45 function evals: %f\n',odecounts)
fprintf('ode113 function evals: %f\n',odecounts1)

figure(1);
subplot(1,3,1);
plot(u(:,1),u(:,2),Color='m',Marker='diamond'); hold on; grid on;
title('Orbital Problem Solution - pc19');
ylabel('y(t)');
xlabel('x(t)');
subplot(1,3,2);
plot(u1(:,1),u1(:,2),Color='c',Marker='diamond');
title('Orbital Problem Solution - ode45');
ylabel('y(t)');
xlabel('x(t)');
subplot(1,3,3);
plot(u2(:,1),u2(:,2),Color='y',Marker='diamond');
title('Orbital Problem Solution - ode113');
ylabel('y(t)');
xlabel('x(t)');

