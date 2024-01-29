%% model aided data fusion
load('Q2.mat')
dt = 0.01;
ax = a(:,1);
ay = a(:,2);
az = a(:,3);
%% measurement 
theta_m = atan2(-ax, sqrt(ay.^2 + az.^2));
theta_sim = [time', theta_m];
%% quad model
A = [0 1;
    0 0];
B = [0; 0.72];
C = [1 0];

Q = 0.05;
R = 1;
%% kalman filter
x = [0; 0];
P = [10 0;
    0 10];
x_hat = zeros(2, length(theta_m));
P_hat = zeros(2, 2, length(theta_m));
K = zeros(2, length(theta_m));
for i = 1:length(theta_m)
    x_hat(:,i) = x;
    P_hat(:,:,i) = P;
    K(:,i) = P*C'/(C*P*C' + R);
    x_dot = A*x + B*u(i);
    x = x + x_dot*dt;
    P = (eye(2) - K(:,i)*C)*P*(eye(2) - K(:,i)*C)' + K(:,i)*R*K(:,i)';
end

%% plot
figure(1)
plot(time, theta_m, 'r', time, x_hat(1,:), 'b')
legend('measurement', 'estimation')
xlabel('time')
ylabel('angle')
title('Kalman filter')
figure(2)
plot(time, x_hat(2,:))
xlabel('time')
ylabel('angular velocity')
title('Kalman filter')

