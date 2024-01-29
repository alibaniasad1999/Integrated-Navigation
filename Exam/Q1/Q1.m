clear;
clc;
%% euler linirization and use ekf for estimation
syms phi theta psi p q r b real
%% add biad to measurement
p_m = p + b;
q_m = q + b;
r_m = r + b;
%% constant
tau = 382;
dt = 0.01;
%% equation
phi_dot = p_m + tan(theta)*(q_m*sin(phi) + r_m*cos(phi))*tan(theta);
theta_dot = q_m*cos(phi) - r_m*sin(phi);
psi_dot = (q_m*sin(phi) + r_m*cos(phi))/cos(theta);
b_dot = -1/tau*b;
%% function
F = [phi_dot;theta_dot;psi_dot;b_dot];
x = [phi;theta;psi;b];
u = [p;q;r];
%% jacobian
A = jacobian(F,x);
B = jacobian(F,u);
%% mkae function
matlabFunction(F,'File','F_func','Vars',{x,u});
matlabFunction(A,'File','A_func','Vars',{x,u});
matlabFunction(B,'File','B_func','Vars',{x,u});
%% linearization
A_num = subs(A,[phi,theta,psi,p,q,r,b],[0,0,0,0,0,0,0]);
B_num = subs(B,[phi,theta,psi,p,q,r,b],[0,0,0,0,0,0,0]);
C = eye(2, 4);
%% observability
O = obsv(A_num,C);
if length(A) == rank(O)
    disp('observable');
else
    disp('not observable');
end 

%% Part b
%% load data
load('Q1.mat');
%% initial
x = [0.0052, 0.0052, 0.0017, 0]';
P = 1e-6*eye(4,4);
Q = diag([0.001,0.001,0.001]);
R = diag([0.005,0.005]);
phi_hat = zeros(length(time),1);
theta_hat = zeros(length(time),1);
psi_hat = zeros(length(time),1);
b_hat = zeros(length(time),1);
u = Gyroscope';
ax = Accelrometer(:, 1);
ay = Accelrometer(:, 2);
az = Accelrometer(:, 3);
phi_m = atan(ay/az);
theta_m = atan2(-ax, sqrt(ay.^2 + az.^2));
%% ekf
for i = 1:length(time)
    %% prediction
    F = F_func(x,u(:,i));
    A = A_func(x,u(:,i));
    B = B_func(x,u(:,i));
    F_matrix = eye(4) + A*dt;
    Gamma = B*dt;
    x = x + F*0.01;
    P = F_matrix*P*F_matrix' + Gamma*Q*Gamma';
    %% update
    y = [phi_m(i);theta_m(i)] - C*x;
    S = C*P*C' + R;
    K = P*C'/S;
    x = x + K*y;
    P = (eye(4,4) - K*C)*P;
    %% save
    phi_hat(i) = x(1);
    theta_hat(i) = x(2);
    psi_hat(i) = x(3);
    b_hat(i) = x(4);
end

%% plot
figure(1);
plot(time, phi_m, 'b');
hold on;
plot(time, phi_hat, 'r');
plot(time, Euler_true(:, 1));
xlabel('time');
ylabel('phi');
legend('measurement', 'estimation', 'true');
% save fig
saveas(gcf,'Q1_phi.png');

figure(2);
plot(time, theta_m, 'b');
hold on;
plot(time, theta_hat, 'r');
plot(time, Euler_true(:, 2));
xlabel('time');
ylabel('theta');
legend('measurement', 'estimation', 'true');
saveas(gcf,'Q1_theta.png');

figure(3);
plot(time, psi_hat, 'r');
hold on;
plot(time, Euler_true(:, 3));
xlabel('time');
ylabel('psi');
legend('estimation', 'true');
saveas(gcf,'Q1_psi.png');

figure(4);
plot(time, b_hat, 'r');
xlabel('time');
ylabel('bias');
saveas(gcf,'Q1_bias.png');
