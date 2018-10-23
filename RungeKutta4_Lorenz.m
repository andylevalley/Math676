clc;                                             
clear all;

% time step and final time
dt = [0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1];
tf = 100;

format long
                
for k = 1:length(dt)
    
    t = 0:dt(k):tf;
    x = zeros(3,length(t)); 
    
    x(1:3,1) = [1;1;1];   
    
    for i=1:(length(t)-1)                              
        k1 = Lorenz(t(i),x(:,i));
        k2 = Lorenz(t(i)+0.5*dt(k),x(:,i)+0.5*dt(k)*k1);
        k3 = Lorenz((t(i)+0.5*dt(k)),(x(:,i)+0.5*dt(k)*k2));
        k4 = Lorenz((t(i)+dt(k)),(x(:,i)+k3*dt(k)));

        x(:,i+1) = x(:,i) + (1/6)*(k1+2*k2+2*k3+k4)*dt(k); 
    end
    
    Error(1:3,k) = x(:,2/dt(k));
    
end

figure(1)
plot3(x(1,:),x(2,:),x(3,:),'k')
title('Lorenz System, $\Delta t = 0.0001$, $\mathbf{X}_0 = [1,1,1]$')
grid on
xlabel('x');
ylabel('y');
zlabel('z');


for i = 1:length(dt)-1
    absError(i) = norm(abs(Error(1:3,i)-Error(1:3,i+1)));
end

dt = [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1];
data = [0,2,4,8];

figure(2)
loglog(dt,absError,'-ok','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
loglog(dt,dt,'-r','MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
loglog(dt,absError(end)/.1^3*dt.^2,'-g','MarkerFaceColor','g','MarkerEdgeColor','g')
hold on
loglog(dt,absError(end)*10^12*dt.^3,'-b','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
loglog(dt,absError(end)*10^(16)*dt.^4,'-m','MarkerFaceColor','m','MarkerEdgeColor','m')
hold on
loglog(dt,absError(end)*10^(20)*dt.^5,'-c','MarkerFaceColor','c','MarkerEdgeColor','c')
axis tight
grid on
title('Cauchy Error');
xlabel('Log($k$)')
ylabel('Log($C_n$)')
legend({'Actual Data','$k$','$k^2$','$k^3$','$k^4$','$k^5$'})



function dxdt = Lorenz(t,x)

sigma = 9;
beta = 1;
rho = 26;

dxdt = [sigma*(x(2)-x(1))
        x(1)*(rho-x(3))-x(2)
        x(1)*x(2)-beta*x(3)];
    
end