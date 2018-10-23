clc
clear all


%% Implement Crank-Nicholoson

Lx = 3/2;
Lt = 4;
dx = .1;
dt = .1;

x = 0:dt:Lx;
Nx = length(x);
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
        

for n = 2:Nt
    b = A*u(1:Nx,n-1);
    u(1:Nx,n) = b\A;
end

figure(1)
surf(x,t,u)
figure(2)
plot(x,u)

%%
Lx = 3/2;
Lt = 4;
dx = 0.01;
dt = 0.01;

t = linspace(0,Lt,100);
x = linspace(0,Lx,100);
n = 1;
f = exp((-4*n*pi^2/9)*t);

plot(x,f)



    
    
