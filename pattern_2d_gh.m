% Code prepared by Kalyan Manna
% For 2D Pattern Formation in Protest Model
% Last updated on 26/04/25

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
g=10;% gamma
d1=0.1;
d2=1;

% Grids
x_min=0;x_max=10;
y_min=0;y_max=10;
t_min=0;t_max=100;
dx=0.05;dy=0.05;dt=0.0001;
x=x_min:dx:x_max;
y=y_min:dy:y_max;
t=t_min:dt:t_max;

Nx=length(x);
Ny=length(y);
Nt=length(t);

% Coexistence steady state
Us=28.7164; Vs=57.4328;

% Initial population distribution
for j=1:Nx
    for k=1:Ny
        u(j,k)=Us+0.001*randn;
        v(j,k)=Vs+0.001*randn;
    end
end

U(1)=mean(mean(u));
V(1)=mean(mean(v));

% Applying FTCS method
for i=2:Nt
    for j=2:Nx-1
        for k=2:Ny-1
            uu(j,k)=u(j,k)+dt*d1*(u(j+1,k)+u(j-1,k)+u(j,k+1)+u(j,k-1)-4*u(j,k))/(dx)^2 ...
                +dt*u(j,k)*(e-v(j,k)+a*u(j,k)/((h^2) +(u(j,k)^2)));
            vv(j,k)=v(j,k)+dt*d2*(v(j+1,k)+v(j-1,k)+v(j,k+1)+v(j,k-1)-4*v(j,k))/(dx)^2 ...
                +dt*(b*u(j,k)-g*v(j,k));
        end
    end
    for j=2:Nx-1
        uu(j,1)=u(j,1)+dt*d1*(u(j+1,1)+u(j-1,1)+u(j,2)-3*u(j,1))/(dx)^2 ...
            +dt*u(j,1)*(e-v(j,1)+a*u(j,1)/((h^2) +(u(j,1)^2)));
        vv(j,1)=v(j,1)+dt*d2*(v(j+1,1)+v(j-1,1)+v(j,2)-3*v(j,1))/(dx)^2 ...
            +dt*(b*u(j,1)-g*v(j,1));
        uu(j,Ny)=u(j,Ny)+dt*d1*(u(j+1,Ny)+u(j-1,Ny)+u(j,Ny-1)-3*u(j,Ny))/(dx)^2 ...
            +dt*u(j,Ny)*(e-v(j,Ny)+a*u(j,Ny)/((h^2) +(u(j,Ny)^2)));
        vv(j,Ny)=v(j,Ny)+dt*d2*(v(j+1,Ny)+v(j-1,Ny)+v(j,Ny-1)-3*v(j,Ny))/(dx)^2 ...
            +dt*(b*u(j,Ny)-g*v(j,Ny));
    end
    for k=2:Ny-1
        uu(1,k)=u(1,k)+dt*d1*(u(2,k)+u(1,k+1)+u(1,k-1)-3*u(1,k))/(dx)^2 ...
            +dt*u(1,k)*(e-v(1,k)+a*u(1,k)/((h^2) +(u(1,k)^2)));
        vv(1,k)=v(1,k)+dt*d2*(v(2,k)+v(1,k+1)+v(1,k-1)-3*v(1,k))/(dx)^2 ...
            +dt*(b*u(1,k)-g*v(1,k));
        uu(Nx,k)=u(Nx,k)+dt*d1*(u(Nx-1,k)+u(Nx,k+1)+u(Nx,k-1)-3*u(Nx,k))/(dx)^2 ...
            +dt*u(Nx,k)*(e-v(Nx,k)+a*u(Nx,k)/((h^2) +(u(Nx,k)^2)));
        vv(Nx,k)=v(Nx,k)+dt*d2*(v(Nx-1,k)+v(Nx,k+1)+v(Nx,k-1)-3*v(Nx,k))/(dx)^2 ...
            +dt*(b*u(Nx,k)-g*v(Nx,k));
    end
    uu(1,1)=u(1,1)+dt*d1*(u(2,1)+u(1,2)-2*u(1,1))/(dx)^2 ...
        +dt*u(1,1)*(e-v(1,1)+a*u(1,1)/((h^2) +(u(1,1)^2)));
    vv(1,1)=v(1,1)+dt*d2*(v(2,1)+v(1,2)-2*v(1,1))/(dx)^2 ...
        +dt*(b*u(1,1)-g*v(1,1));
    uu(1,Ny)=u(1,Ny)+dt*d1*(u(2,Ny)+u(1,Ny-1)-2*u(1,Ny))/(dx)^2 ...
        +dt*u(1,Ny)*(e-v(1,Ny)+a*u(1,Ny)/((h^2) +(u(1,Ny)^2)));
    vv(1,Ny)=v(1,Ny)+dt*d2*(v(2,Ny)+v(1,Ny-1)-2*v(1,Ny))/(dx)^2 ...
        +dt*(b*u(1,Ny)-g*v(1,Ny));
    uu(Nx,1)=u(Nx,1)+dt*d1*(u(Nx-1,1)+u(Nx,2)-2*u(Nx,1))/(dx)^2 ...
        +dt*u(Nx,1)*(e-v(Nx,1)+a*u(Nx,1)/((h^2) +(u(Nx,1)^2)));
    vv(Nx,1)=v(Nx,1)+dt*d2*(v(Nx-1,1)+v(Nx,2)-2*v(Nx,1))/(dx)^2 ...
        +dt*(b*u(Nx,1)-g*v(Nx,1));
    uu(Nx,Ny)=u(Nx,Ny)+dt*d1*(u(Nx-1,Ny)+u(Nx,Ny-1)-2*u(Nx,Ny))/(dx)^2 ...
        +dt*u(Nx,Ny)*(e-v(Nx,Ny)+a*u(Nx,Ny)/((h^2) +(u(Nx,Ny)^2)));
    vv(Nx,Ny)=v(Nx,Ny)+dt*d2*(v(Nx-1,Ny)+v(Nx,Ny-1)-2*v(Nx,Ny))/(dx)^2 ...
        +dt*(b*u(Nx,Ny)-g*v(Nx,Ny));
    
    u=uu;
    v=vv;
    
    U(i)=mean(mean(u));
    V(i)=mean(mean(v));
    
    %Figure
    if mod(i,100)==0
        figure(1)
        pcolor(x,y,u); shading interp; colorbar;
        pause(0.00001)
        figure(2)
        pcolor(x,y,v); shading interp; colorbar;
        pause(0.00001)
        i=i
    end
end