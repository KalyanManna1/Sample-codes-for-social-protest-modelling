% Code prepared by Kalyan Manna
% For 1D Pattern Formation in Protest Model
% Last updated on 05/03/25

clear all;
clc;
clf;
close all;

% Format
format long e

set(0,'DefaultAxesFontSize',20);

% Parameter values
e=1;% epsilon
a=20000;
h=110;
b=20;% beta
g=13;% gamma
d1=0.2;
d2=1;

% Grids
x_min=0;x_max=10;
t_min=0;t_max=2000;
dx=0.05;dt=0.0001;
x=x_min:dx:x_max;
t=t_min:dt:t_max;

Nx=length(x);
Nt=length(t);
mid=(Nx-1)/2+1;

% Coexistence steady state
Us=33.9652; Vs=52.2542;

% Initial population distribution
for j=1:Nx
    u(j)=Us+0.001*randn;
    v(j)=Vs+0.001*randn;
end

Uav(1)=mean(u);
Vav(1)=mean(v);

UU=[];
UUU=[];
VV=[];

for i=2:Nt
    for j=2:Nx-1
        uu(j)=u(j)+d1*(dt/(dx)^2)*(u(j-1)+u(j+1)-2*u(j))+...
            dt*u(j)*(e-v(j)+a*u(j)/((h^2)+(u(j)^2)));
        vv(j)=v(j)+d2*(dt/(dx)^2)*(v(j-1)+v(j+1)-2*v(j))+...
            dt*(b*u(j)-g*v(j));
    end
        
        uu(1)=u(1)+d1*(dt/(dx)^2)*(u(2)-u(1))+...
            dt*u(1)*(e-v(1)+a*u(1)/((h^2)+(u(1)^2)));
        uu(Nx)=u(Nx)+d1*(dt/(dx)^2)*(u(Nx-1)-u(Nx))+...
            dt*u(Nx)*(e-v(Nx)+a*u(Nx)/((h^2)+(u(Nx)^2)));
        vv(1)=v(1)+d2*(dt/(dx)^2)*(v(2)-v(1))+...
            dt*(b*u(1)-g*v(1));
        vv(Nx)=v(Nx)+d2*(dt/(dx)^2)*(v(Nx-1)-v(Nx))+...
            dt*(b*u(Nx)-g*v(Nx));
        
        u=uu;
        v=vv;
        
        Uav(i)=mean(u);
        Vav(i)=mean(v);

        % Figure
        if mod(i,1000)==0
            figure(1)
            plot(x,u);
            pause(0.00001)
            UU=[UU; u];
            VV=[VV; v];
%             if i>=900000
%                 UUU=[UUU; u];
%             end
            i=i
        end
end

figure(2)
plot(Uav,Vav);
hold on
plot(Us,Vs,'r.','Markersize',20);
xlabel('$\langle U\rangle$','Interpreter','Latex');
ylabel('$\langle V\rangle$','Interpreter','Latex');
% set(gca, 'xlim', [32 40]);
% set(gca, 'ylim', [50 60]);
hold off

figure(3)
pcolor(UU); shading interp; 
colorbar;
%set(cbh,'YTick',[2:0.1:3])
%c.Limits = [0 3];
%c.Ticks = [2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3];
ax = gca;
ax.XTick=[];
ax.YTick=[];
xlabel('$x \rightarrow$','Interpreter','Latex');
ylabel('$t \rightarrow$','Interpreter','Latex');

figure(4)
fig = figure;
left_color = [0 0 1];
right_color = [1 0 1];
% right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(x,u,'linewidth',2);
xlabel('$x \rightarrow$','Interpreter','Latex');
ylabel('$U$','Interpreter','Latex');
set(gca, 'xlim', [0 10]);
hold on
yyaxis right
plot(x,v,'linewidth',2);
ylabel('$V$','Interpreter','Latex');
hold off