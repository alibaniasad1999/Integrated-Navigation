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
%% ekf %%
%% constant %%
Q = 1*eye(5);
R = 0.1;
%% initial %%
x0 = [0;0;0;0;0];
P0 = 0.1*eye(5);
%% loas data u and z %%
load('u.mat');
load('z.mat');
dt = 0.01;
T = 0:dt:10;
x_hat_save = zeros(5,length(T));
P_save = zeros(5,5,length(T));
%% ekf %%
x_hat = x0;
P = P0;
for i = 1:length(T)
    % prediction
    x_hat = x_hat + dt*subs(f,x,x_hat);
    A = jacobian(f,x);
    P = P + dt*(A*P + P*A' + Q);
    % update
    y = z(i) - x_hat(1);
    C = jacobian(y,x);
    S = C*P*C' + R;
    K = P*C'/S;
    x_hat = x_hat + K*y;
    P = P - K*S*K';
    % save
    x_hat_save(:,i) = x_hat;
    P_save(:,:,i) = P;
end
%% plot %%
figure(1)
plot(T,x_hat_save(1,:),'r','linewidth',2);
hold on;
plot(T,z,'b','linewidth',2);
xlabel('time(s)');
ylabel('height(m)');
legend('estimate','measurement');
figure(2)
plot(T,x_hat_save(2,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('velocity(m/s)');
figure(3)
plot(T,x_hat_save(3,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('bias');
figure(4)
plot(T,x_hat_save(4,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('a0');
figure(5)
plot(T,x_hat_save(5,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('a1');
figure(6)
plot(T,P_save(1,1,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('P11');
figure(7)
plot(T,P_save(2,2,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('P22');
figure(8)
plot(T,P_save(3,3,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('P33');
figure(9)
plot(T,P_save(4,4,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('P44');
figure(10)
plot(T,P_save(5,5,:),'r','linewidth',2);
xlabel('time(s)');
ylabel('P55');

