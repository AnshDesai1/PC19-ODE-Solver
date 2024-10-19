clear all; close all; clc;
y0=cos(1)-1;
t0=0;
tf=6;
tol=1e-5;
tspan=[t0 tf];
opt=odeset('AbsTol',.8e-8,'RelTol',1.01e-8);

u=@(t)exp(sin(12*t)).*(cos((t-1).^2)-1);
u1=@(t)12*exp(sin(12*t)).*(-1+cos((-1+t).^2)).*cos(12*t)-2*exp(sin(12*t)).*(-1+t).*sin((-1+t).^2);
f=@(t,y)y.^2+u1(t)-u(t).^2;

[y,t,counter]=pc19(f,t0,y0,tf,tol,9);
SOL=ode45(f,tspan,y0,opt);
SOL1=ode113(f,tspan,y0,opt);
t1=SOL.x';
y1=SOL.y';
t2=SOL.x';
y2=SOL.y';
odecounts=SOL.stats.nfevals;
odecounts1=SOL1.stats.nfevals;

fprintf('pc19 function evals: %f\n',counter)
fprintf('ode45 function evals: %f\n',odecounts)
fprintf('ode113 function evals: %f\n',odecounts1)

figure(1);
subplot(1,2,1);

plot(t,abs(y-u(t)),'c', t1,abs(y1-u(t1)),'r',t2,abs(y2-u(t2)),'g');
title('Wiggly Problem Error');
ylabel('|y-u(t)|');
xlabel('t');
legend('pc19','ode45','ode113')

subplot(1,2,2);
plot(t,y,Color='m');
title('Wiggly Problem Solution');
ylabel('Solution y(t)');
xlabel('t');