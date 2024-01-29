clc;
clear;
lin;
%% INS equation with the source term
%% load data
load('Q3.mat');
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
position = in_profile_Master(1, 2:4)';
velocity = in_profile_Master(1, 5:7)';
attitude = in_profile_Master(1, 8:10)';
position_rate = zeros(3, 1);
velocity_rate = zeros(3, 1);
attitude_rate = zeros(3, 1);
position_array = zeros(length(Accelerometer_Slave), 3);
velocity_array = zeros(length(Accelerometer_Slave), 3);
attitude_array = zeros(length(Accelerometer_Slave), 3);
P = 10*eye(9);
Q = 0.1*eye(6);
R = 0.1*eye(3);
x = [position; velocity; attitude];

%% load master data
out = sim('INS_GPS_sim');
Master = out.simout.Data;
Cbn = Cbn_calculation(attitude);
for i = 1:length(Accelerometer_Slave)
    R_m = a*(1-ee)/((1-ee*(sin(position(1)))^2)^(3/2));
    R_p = a/((1-ee*(sin(position(1)))^2)^(1/2));
    R_G = sqrt(R_m*R_p);
    %% Calculating Force at body frame
    Fb = Accelerometer_Slave(i, 2:4)';
    %% Calculating Angular rate at body frame
    wib = Gyroscope_Slave(i, 2:4)';
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
    %% save data
    position_array(i, :) = position';
    velocity_array(i, :) = velocity';
    %% Kalman filter
    V_SN = Master(i, 5:7);
    % using the velocity from the master as the measurement
    velocity_master = V_SN';
    % x = [position; velocity; attitude];
    % % delete bias form A matrix
    % A = A_num(1:9, 1:9);
    % F = eye(9) + A*dt;
    % B = B_num(1:9, 1:6);
    % Gamma = B*dt;
    % P = F*P*F' + Gamma*Q*Gamma';
    % C = [zeros(3, 3), eye(3), zeros(3, 3)];
    % K = P*C'/(C*P*C'+R);
    % x = x + K*(velocity_master - C*x);
    % P = (eye(9) - K*C)*P;
    % position = x(1:3);
    % velocity = x(4:6);
    % attitude = x(7:9);
    % attitude_array(i, 2) = asin(-Cbn(3, 1)); % theta
    % attitude_array(i, 1) = acos(Cbn(3, 3)/cos(attitude_array(i, 2)))*sign(Cbn(3, 2)); % phi
    % attitude_array(i, 3) = acos(Cbn(1, 1)/cos(attitude_array(i, 2)))*sign(Cbn(2, 1)); % psi
    % DCM to Euler angle matlab function
    attitude_array(i, :) = attitude';
end
%%
figure(1)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), position_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('latitude', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(2)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), position_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('longitude', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(3)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), position_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('altitude', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(4)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), velocity_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in x', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(5)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), velocity_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in y', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(6)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), velocity_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('velocity in z', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(7)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), attitude_array(:, 1), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\phi$', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(8)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), attitude_array(:, 2), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\theta$', 'interpreter', 'latex', 'FontSize', 24);
title('');

figure(9)
set(gca, 'FontSize', 16)
hold on;
plot(Accelerometer_Slave(1:end, 1), attitude_array(:, 3), 'LineWidth', 2,...
    'Color','k', 'linestyle', '--');
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman');
xlabel('time($\sec)$', 'interpreter', 'latex', 'FontSize', 24);
ylabel('$\psi$', 'interpreter', 'latex', 'FontSize', 24);
title('');







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
