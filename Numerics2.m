clc
clear all


%% Implement Crank-Nicholoson

Lx = 3/2;
Lt = 4;
Nx = 200;
Nt = 200;

x = linspace(0,Lx,Nx);
dx = x(2) - x(1);
t = linspace(0,Lt,Nx);
Nt = length(t);
dt = t(2) - t(1);

lambda = dt/dx^2;

A = spdiags(repmat([-(1/2)*lambda, 1+2+lambda, -(1/2)*lambda],[Nx 1]),-1:1,Nx,Nx);
A(1,1) = 1;
A(1,2) = 0;
A(Nx,Nx-1) = 0;
A(Nx,Nx) = 1;

u = zeros(Nx,Nt);
for i = 1:Nx
    if x(i) < 1/2
        u(i,1) = x(i);
    elseif (x(i) >= 1/2) && (x(i) <= 1)
        u(i,1) = 0;
    elseif x(i) > 1
        u(i,1) = -x(i) + 3/2;
    end
end


plot(x,u(1:Nx,1))
title(sprintf('t = %f',t(1)))
for n = 2:Nt
    b = A*u(1:Nx,n-1);
    u(1:Nx,n) = b\A;
    
    figure(1)
    plot(x,u(1:Nx,n));
    title(sprintf('t = %f',t(n)))
    drawnow;
end

figure(1)
plot(x,u(:,end))

%%
n = 1;
u_exact = u;

for n = 2:Nt
    u_exact(:,n) = u_exact(:,n-1)*exp((-4*n*pi^2/9)*dt);
end

figure(2)
plot(x,u_exact(:,end))



    
    
