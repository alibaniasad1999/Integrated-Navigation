clc;
clear;
%% state space %%
syms x1 x2 x3 x4 x5 u
% states 1 ---> height
% states 2 ---> velocity
% states 3 ---> bias
% states 4 ---> a0
% states 5 ---> a1
% constant
h_r = 0.013;
R = 0.0333;
m = 0.063;
%% state space
f1 = x2;
z = x1 + x3;
h_bar = (-z + h_r)/(2.5*R);
K_G = x4 + x5*(2/h_bar)^2;
f2 = -9.81 - K_G/m*u;
f3 = 0;
f4 = 0;
f5 = 0;
f = [f1;f2;f3;f4;f5];
x = [x1;x2;x3;x4;x5];
%% particle filter %%
% initial condition
x0 = [0;0;0;0;0];
% noise
Q = diag([0.0001,0.0001,0.0001,0.0001,0.0001]);
R = 0.0001;
% time
dt = 0.01;
t = 0:dt:10;
% input
u = 0.5*ones(1,length(t));
% measurement
y = zeros(1,length(t));
% initial condition
x_hat = zeros(5,length(t));
x_hat(:,1) = x0;
P = zeros(5,5,length(t));
P(:,:,1) = diag([0.0001,0.0001,0.0001,0.0001,0.0001]);
% particle
N = 1000;
x_hat_particle = zeros(5,N,length(t));
x_hat_particle(:,1,:) = repmat(x0,1,N);
w = zeros(N,length(t));
w(:,1) = 1/N*ones(N,1);
% particle filter
for i = 2:length(t)
    % state space
    x_hat(:,i) = x_hat(:,i-1) + dt*f;
    % measurement
    y(i) = x_hat(1,i) + sqrt(R)*randn;
    % particle filter
    for j = 1:N
        % state space
        x_hat_particle(:,j,i) = x_hat_particle(:,j,i-1) + dt*f + sqrt(Q)*randn(5,1);
        % measurement
        y_particle = x_hat_particle(1,j,i) + sqrt(R)*randn;
        % weight
        w(j,i) = w(j,i-1)*exp(-1/2*(y(i)-y_particle)^2/R);
    end
    % normalize
    w(:,i) = w(:,i)/sum(w(:,i));
    % resampling
    index = randsample(1:N,N,true,w(:,i));
    x_hat_particle(:,:,i) = x_hat_particle(:,index,i);
    w(:,i) = 1/N*ones(N,1);
end
%% plot
figure(1)
plot(t,x_hat(1,:),'b','linewidth',2)
hold on
plot(t,y,'r','linewidth',2)
plot(t,x_hat_particle(1,:,end),'k','linewidth',2)
xlabel('time')
ylabel('height')
legend('state space','measurement','particle filter')
figure(2)
plot(t,x_hat(2,:),'b','linewidth',2)
hold on
plot(t,x_hat_particle(2,:,end),'k','linewidth',2)
xlabel('time')
ylabel('velocity')
legend('state space','particle filter')
figure(3)
plot(t,x_hat(3,:),'b','linewidth',2)
hold on
plot(t,x_hat_particle(3,:,end),'k','linewidth',2)
xlabel('time')
ylabel('bias')
legend('state space','particle filter')
figure(4)
plot(t,x_hat(4,:),'b','linewidth',2)
hold on
plot(t,x_hat_particle(4,:,end),'k','linewidth',2)
xlabel('time')
ylabel('a0')
legend('state space','particle filter')
figure(5)
plot(t,x_hat(5,:),'b','linewidth',2)
hold on
plot(t,x_hat_particle(5,:,end),'k','linewidth',2)
xlabel('time')
ylabel('a1')
legend('state space','particle filter')
figure(6)
plot(t,w(:,end),'k','linewidth',2)
xlabel('time')
ylabel('weight')
legend('particle filter')

