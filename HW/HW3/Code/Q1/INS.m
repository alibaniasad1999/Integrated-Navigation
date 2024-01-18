clc;
clear;
global P_GPS
P_GPS = 0.01*eye(9);
%% INS equation with the source term
%% load data
load('DataSet_02.mat');
%% parameters:
%% Sample rate
dt = 0.01; %Sample time
%% Earth Parameters
f = 1/298.2572236;
e = 0.0818191908426;
R_0 = 6378137;
ee = e*e;
Omega_e = 7.29211651452621e-5; %Earth rotation rate
a = 6378178;
%% initial conditions
position = in_profile(1, 2:4)';
velocity = in_profile(1, 5:7)';
attitude = in_profile(1, 8:10)';
position_rate = zeros(3, 1);
velocity_rate = zeros(3, 1);
attitude_rate = zeros(3, 1);
position_array = zeros(length(IMU_true), 3);
velocity_array = zeros(length(IMU_true), 3);
attitude_array = zeros(length(IMU_true), 3);
Cbn = Cbn_calculation(attitude);
for i = 1:length(IMU_true)
    R_m = a*(1-ee)/((1-ee*(sin(position(1)))^2)^(3/2));
    R_p = a/((1-ee*(sin(position(1)))^2)^(1/2));
    R_G = sqrt(R_m*R_p);
    %% Calculating Force at body frame
    Fb = IMU_true(i, 2:4)';
    %% Calculating Angular rate at body frame
    wib = IMU_true(i, 5:7)';
    omega_i_e_n = [Omega_e*cos(position(1)); 0;...
        -Omega_e*sin(position(1))];
    omega_e_n_n = [velocity(2)/(R_p+position(3));...
        -velocity(1)/(R_m+position(3));...
        -velocity(2)*tan(position(1))/(R_m+position(3))];
    omega_i_n_n = omega_i_e_n + omega_e_n_n;
    Cbn = Cbn_calculation(attitude);
    omega_n_b_b = wib - Cbn'*omega_i_n_n;
    %% Calculating Angular rate at navigation frame
    attitude_rate = [1, sin(attitude(1))*tan(attitude(2)), cos(attitude(1))*tan(attitude(2));
        0, cos(attitude(1)), -sin(attitude(1));
        0, sin(attitude(1))/cos(attitude(2)), cos(attitude(1))/cos(attitude(2))]*omega_n_b_b;
    %% Calculating Acceleration at body frame
    Fn = Cbn*Fb;
    %% Calculating Acceleration at navigation frame
    Fg = [0; 0; 9.780318*(1+5.3024e-3*(sin(position(1)))^2-0.0000059*(sin(2*position(1))))];
    Fg = Fg / (1 + position(3)/R_G)^2;
    %% correction terms
    correction = cross(2*omega_i_e_n+omega_e_n_n, velocity);
    %% Calculating Acceleration at navigation frame
    F = Fn + Fg - correction;
    %% Calculating Velocity at navigation frame
    velocity_rate = F;
    %% Calculating Position at navigation frame
    Cen = [cos(position(1)), -sin(position(1))*sin(position(2)), sin(position(1))*cos(position(2));
        0, cos(position(2)), sin(position(2));
        -sin(position(1)), -cos(position(1))*sin(position(2)), cos(position(1))*cos(position(2))];
    Cen_rate = - skew_symmetric(omega_e_n_n)*Cen;
    %% integral
    attitude = attitude + attitude_rate*dt;
    velocity = velocity + velocity_rate*dt;
    Cen = Cen + Cen_rate*dt;
    position(1) = acos(Cen(1, 1));
    position(2) = asin(Cen(2, 3));
    position(3) = position(3) - velocity(3)*dt;
    %% add GPS data
    if mod(i,100) == 0
        g = 9.7803267714;
        R_m = 6378137;
        R_p = 6356752.3142;
        R_g = 6371008.7714;
        Omega_e = 7.292116e-5;
        x = [position; velocity; attitude];
        u = [Fb; wib];
        A = A_func(x, u);
        B = B_func(x, u);
        F = eye(9) + dt*A;
        % x = F*[position; velocity; attitude];
        % position = x(1:3);
        % velocity = x(4:6);
        % attitude = x(7:9);
        GPS = GPS_meas(i/100, 2:end)';
        d_GPS = GPS - [position; velocity];
        [d_position, d_velocity, d_attitude] = GPS_INS_KF(position, velocity, attitude, A, B, d_GPS, dt);
        position = position + d_position;
        velocity = velocity + d_velocity;
        attitude = attitude + d_attitude;
    end
    %% save data
    position_array(i, :) = position';
    velocity_array(i, :) = velocity';
    attitude_array(i, :) = attitude';
end
% figure(1)
% set(gca, 'FontSize', 16)
% hold on;
% plot(in_profile(:, 1), in_profile(:, 2), 'LineWidth', 2, 'Color','r');
% plot(in_profile(1:end-1, 1), position_array(:, 1), 'LineWidth', 2,...
%     'Color','k', 'linestyle', '--');
% legend('true', 'estimated');
% set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
% xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
% ylabel('latitude', 'interpreter', 'latex', 'FontSize', 24);
% title('');
% axis tight
% [dir_state, ~, ~] = mkdir('../../../Figure/Q5');
% if dir_state
%     print('../../../Figure/Q5/latitude','-depsc');
% else
%     fprintf("Ooooooops\n")
% end

figure(1)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 2), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), position_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('latitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(2)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 3), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), position_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('longitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(3)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 4), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), position_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('altitude', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(4)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 5), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), velocity_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in x', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(5)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 6), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), velocity_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in y', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(6)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 7), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), velocity_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in z', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(7)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 8), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), attitude_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\phi$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(8)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 9), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), attitude_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time(sec)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\theta$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight

figure(9)
set(gca, 'FontSize', 16)
hold on;
plot(in_profile(:, 1), in_profile(:, 10), 'LineWidth', 2, 'Color','r');
plot(in_profile(1:end-1, 1), attitude_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');        
legend('true', 'estimated');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time(sec)', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\psi$', 'interpreter', 'latex', 'FontSize', 24);
title('');
axis tight




%% Cbn calculationz
function Cbn = Cbn_calculation(attitude)
    phi = attitude(1);
    theta = attitude(2);
    psi = attitude(3);
    Cnb = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
        -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi),...
        cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi),...
        sin(phi)*cos(theta);
        sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)...
        , -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi),...
        cos(phi)*cos(theta)];
    Cbn = Cnb';
end
%% skew symmetric matrix
function S = skew_symmetric(w)
    S = [0, -w(3), w(2);
        w(3), 0, -w(1);
        -w(2), w(1), 0];
end


%% GPS INS kalman filter
function [position, velocity, attitude] = GPS_INS_KF(~, ~, ~, A, B, GPS, dt)
    %% EKF
    global P_GPS
    %% state transition matrix
    F = eye(9) + A*dt;
    %% control matrix
    G = B*dt;
    %% covariance matrix
    Q = 0.01*eye(6);
    %% state prediction
    x = zeros(9,1);   
    %% covariance prediction
    P_GPS = F*P_GPS*F' + G*Q*G';
    %% measurement matrix
    H = [eye(3), zeros(3, 6);
        zeros(3, 3), eye(3), zeros(3, 3)];
    %% measurement prediction
    z = GPS;
    %% measurement covariance
    R = 1*eye(6);
    %% Kalman gain
    K = P_GPS*H'/(H*P_GPS*H'+R);
    %% state update
    x = x + K*(z);
    %% covariance update
    P_GPS = (eye(9) - K*H)*P_GPS;
    %% output
    position = x(1:3);
    velocity = x(4:6);
    attitude = x(7:9);
end



