clc;
clear;
%% state space %%
syms x1 x2 x3 x4 x5 u real % x1,x2,x3,x4,x5 are real
% states 1 ---> height
% states 2 ---> velocity
% states 3 ---> bias
% states 4 ---> a0
% states 5 ---> a1
% constant
h_r = 0.013;
R_r = 0.0333;
m = 0.063;
%% state space
f1 = x2;
zm = x1 + x3;
h_bar = (-zm + h_r)/(2.5*R_r);
K_G = x4 + x5*(2/h_bar)^2;
f2 = 9.81 - K_G/m*u;
f3 = 0;
f4 = 0;
f5 = 0;
f = [f1;f2;f3;f4;f5];
% make function for f
x = [x1;x2;x3;x4;x5];
matlabFunction(f,'File','f_func','Vars',{x,u});
%% ekf %%
%% initial %%
x0 = [-2;0;0.2;0.9621;0.03794];
P0 = 1e-5*eye(5);
out = sim('Quad_IGE_Landing');
%% constant %%
Q = 1e-9*eye(5);
R = 1e-1;
%% load data u and z %%
%% run simulink %%
u_in = out.u.Data;
zm = out.zm.Data; % measurement
dt = 0.01;
T = 0:dt:20;
x_hat_save = zeros(5,length(T));
P_save = zeros(5,5,length(T));
%% ekf %%
x_hat = x0;
P = P0;
A = jacobian(f,x);
% make A function
matlabFunction(A,'File','A_func','Vars',{x,u});
A_s = A_func(x_hat,u_in(1));
F = eye(5) + dt*A_s;
h_hat_array = zeros(5,length(T));
for i = 1:length(T)
    % prediction
    x_dot = f_func(x_hat,u_in(1));
    x_hat = x_hat + dt*x_dot;
    % save h_hat 
    h_hat_array(:,i) = double(vpa(x_hat));
    P = F*P*F' + Q;
    % update
    H = [1 0 0 0 0];
    z = H*x_hat;
    K = P*H'/(H*P*H' + R);
    x_hat = x_hat + K*(zm(i) - z);
    P = (eye(5) - K*H)*P*(eye(5) - K*H)' + K*R*K';
    % save
    x_hat_save(:,i) = double(vpa(x_hat));
    P_save(:,:,i) = double(vpa(P));
    % write bar progress of simulation to command window
    clc;
    fprintf('%.2f%%\n',i/length(T)*100);
end
%% plot
figure(1)
plot(T,x_hat_save(1,:),'LineWidth',2)
hold on
plot(T,out.zm.Data,'LineWidth',2)
plot(T,h_hat_array(1,:),'LineWidth',2)
plot(T,out.z.Data,'LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Height [m]')
legend('Estimation','True')
figure(2)
plot(T,x_hat_save(2,:),'LineWidth',2)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
figure(3)
plot(T,x_hat_save(3,:),'LineWidth',2)
xlabel('Time [s]')
ylabel('Bias [m]')
grid on
figure(4)
plot(T,x_hat_save(4,:),'LineWidth',2)
xlabel('Time [s]')
ylabel('a0 [m/s^2]')
grid on
figure(5)
plot(T,x_hat_save(5,:),'LineWidth',2)
xlabel('Time [s]')
ylabel('a1 [m/s^2]')
grid on
%% plot P
figure(6)
subplot(2,3,1)
plot(T,squeeze(P_save(1,1,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P11')
grid on
subplot(2,3,2)
plot(T,squeeze(P_save(2,2,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P22')
grid on
subplot(2,3,3)
plot(T,squeeze(P_save(3,3,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P33')
grid on
subplot(2,3,4)
plot(T,squeeze(P_save(4,4,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P44')
grid on
subplot(2,3,5)
plot(T,squeeze(P_save(5,5,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P55')
grid on
subplot(2,3,6)
plot(T,squeeze(P_save(1,2,:)),'LineWidth',2)
xlabel('Time [s]')
ylabel('P12')
grid on

