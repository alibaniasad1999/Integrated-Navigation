clc;
clear;
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
    %% save data
    position_array(i, :) = position';
    velocity_array(i, :) = velocity';
    % attitude_array(i, 2) = asin(-Cbn(3, 1)); % theta
    % attitude_array(i, 1) = acos(Cbn(3, 3)/cos(attitude_array(i, 2)))*sign(Cbn(3, 2)); % phi
    % attitude_array(i, 3) = acos(Cbn(1, 1)/cos(attitude_array(i, 2)))*sign(Cbn(2, 1)); % psi
    % DCM to Euler angle matlab function
    attitude_array(i, :) = flip(rotm2eul(Cbn, 'ZYX'));
end

% plot phi
figure(1)
plot(attitude_array(:, 1))
hold on
plot(in_profile(:, 8))
title('phi')
legend('estimated', 'true')
% plot theta
figure(2)
plot(attitude_array(:, 2))
hold on
plot(in_profile(:, 9))
title('theta')
legend('estimated', 'true')
% plot psi
figure(3)
plot(attitude_array(:, 3))
hold on
plot(in_profile(:, 10))
title('psi')
legend('estimated', 'true')
% plot position
figure(4)
plot(position_array(:, 1))
hold on
plot(in_profile(:, 2))
title('latitude')
legend('estimated', 'true')
figure(5)
plot(position_array(:, 2))
hold on
plot(in_profile(:, 3))
title('longitude')
legend('estimated', 'true')
figure(6)
plot(position_array(:, 3))
hold on
plot(in_profile(:, 4))
title('altitude')
legend('estimated', 'true')
% plot velocity
figure(7)
plot(velocity_array(:, 1))
hold on
plot(in_profile(:, 5))
title('velocity x')
legend('estimated', 'true')
figure(8)
plot(velocity_array(:, 2))
hold on
plot(in_profile(:, 6))
title('velocity y')
legend('estimated', 'true')
figure(9)
plot(velocity_array(:, 3))
hold on
plot(in_profile(:, 7))
title('velocity z')
legend('estimated', 'true')

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
