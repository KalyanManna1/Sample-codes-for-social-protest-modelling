% Code prepared by Kalyan Manna 
% For Temporal bifurcation diagram for protest model
% Last updated on 20/02/2025

function protest_phase_1

clear all;
clc;
clf;
close all;

% Global parameters
global e a h b g

% Format
format long e

set(0,'DefaultAxesFontSize',20);

% Parameter values
e=1;% epsilon
a=20000;
h=110;
b=20;% beta
g=18;% gamma

% Time range
no_of_days=5000;
t0=0;
tf=no_of_days;
dt=0.01;
rt=tf:-0.01:t0;

% Initial conditions
init1=[0 150];
init2=[0.0001 0];
init3=[78 86];

%ode45 solution to the system using the three different initial conditions
options=odeset('RelTol',1e-4,'AbsTol',1e-5);
[T_ODE451,Y_ODE451]=ode45(@ode45_eqns,[t0:dt:tf],init1,options);
%[T_ODE452,Y_ODE452]=ode45(@ode45_eqns,[t0:0.001:10],init2,options);
[T_ODE453,Y_ODE453]=ode45(@ode45_eqns,[t0:0.0001:500],init3,options);
K1=Y_ODE453(:,1);
K2=Y_ODE453(:,2);
init4=[K1(end) K2(end)];
[T_ODE454,Y_ODE454]=ode45(@ode45_eqns,[t0:0.0001:50],init4,options);

% Figure
plot(Y_ODE451(:,1),Y_ODE451(:,2),'color',[0.5 0.5 0.5],'linewidth',2);
hold on
% plot(Y_ODE452(:,1),Y_ODE452(:,2),'color',[0.5 0.5 0.5],'linewidth',2);
% hold on
plot(Y_ODE453(:,1),Y_ODE453(:,2),'color','b','linewidth',2);
hold on
plot(Y_ODE454(:,1),Y_ODE454(:,2),'color','g','linewidth',2);
hold on
plot(0,0,'r.','Markersize',30);
hold on
plot(78.1644,86.8493,'r.','Markersize',30);
xlabel('$U$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca, 'xlim', [-10 200]);
set(gca, 'ylim', [-10 150]);
grid on
text(5,-5,'$E_{0}$','Interpreter','Latex','Color','[0.4940 0.1840 0.5560]','FontSize',14)
text(78,72,'$E_{*}$','Interpreter','Latex','Color','[0.4940 0.1840 0.5560]','FontSize',14)

% Functions for ode45
function dy = ode45_eqns(t,y)

global e a h b g

dy=zeros(2,1);

dy(1)=e*y(1)-y(1)*y(2)+a*y(1)*y(1)/(h^2 +y(1)^2);
dy(2)=b*y(1)-g*y(2);

%End of code