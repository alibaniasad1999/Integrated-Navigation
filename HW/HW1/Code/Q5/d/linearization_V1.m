clc;
clear;
%% linearization of INS equations using MATLAB symbolic toolbox
syms phi theta psi v_n v_e v_d L lambda h real
syms w_x w_y w_z a_x a_y a_z real
% constants
f = 1/298.2572236;
e = 0.0818191908426;
dt = 0.01;
R_0 = 6378137;
ee = e*e;
Omega_e = 7.29211651452621e-5; %Earth rotation rate
a = 6378178;

% state vector
x = [L lambda h v_n v_e v_d phi theta psi]';
% input vector
u = [a_x a_y a_z w_x w_y w_z]';
% state vector derivative
position = [L, lambda, h];
velocity = [v_n, v_e, v_d];
attitude = [phi, theta, psi];
R_m = a*(1-ee)/((1-ee*(sin(position(1)))^2)^(3/2));
R_p = a/((1-ee*(sin(position(1)))^2)^(1/2));
R_G = sqrt(R_m*R_p);
Fb = [a_x, a_y, a_z]';
wib = [w_x, w_y, w_z]';
omega_i_e_n = [Omega_e*cos(position(1)); 0;...
    -Omega_e*sin(position(1))];
omega_e_n_n = [velocity(2)/(R_p+position(3));...
    -velocity(1)/(R_m+position(3));...
    -velocity(2)*tan(position(1))/(R_m+position(3))];
omega_i_n_n = omega_i_e_n + omega_e_n_n;
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
F = Fn + Fg - correction';
%% Calculating Velocity at navigation frame
velocity_rate = F;
position_rate = [velocity(1)/(R_m+position(3));...
    velocity(2)/((R_p+position(3))*cos(position(1)));
    -velocity(3)];
%% Calculating state vector derivative
x_dot = [position_rate; velocity_rate; attitude_rate];
%% Calculating Jacobian matrix
A = jacobian(x_dot, x);
B = jacobian(x_dot, u);
%% save A and B to function
matlabFunction(simplify(A, 'step',100), 'File', 'A_func', 'Vars', {x, u});
matlabFunction(simplify(B, 'step',100), 'File', 'B_func', 'Vars', {x, u});
%% load data
load('DataSet_02.mat')
position_num = in_profile(1, 2:4);
velocity_num = in_profile(1, 5:7);
attitude_num = in_profile(1, 8:10);
%% initial state
x_num = [position_num, velocity_num, attitude_num]';
u_num = IMU_true(1, 2:7)';
position_array = zeros(length(IMU_true), 3);
velocity_array = zeros(length(IMU_true), 3);
attitude_array = zeros(length(IMU_true), 3);
out_profile = zeros(length(IMU_true), 9);
%% run simulation
for i = 1:length(IMU_true)
    u_num = IMU_true(i, 2:7)';
    A_num = A_func(x_num, u_num);
    B_num = B_func(x_num, u_num);
    x_dot_num = A_num*x_num + B_num*u_num;
    x_num = x_num + x_dot_num*dt;
    out_profile(i, :) = x_num';
end
%% plot
figure(1);
plot(IMU_true(:, 1), out_profile(:, 7), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 8), 'b');
hold off;

figure(2);
plot(IMU_true(:, 1), out_profile(:, 8), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 9), 'b');
hold off;

figure(3);
plot(IMU_true(:, 1), out_profile(:, 9), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 10), 'b');
hold off;

figure(4);
plot(IMU_true(:, 1), out_profile(:, 4), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 5), 'b');
hold off;

figure(5);
plot(IMU_true(:, 1), out_profile(:, 5), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 6), 'b');
hold off;

figure(6);
plot(IMU_true(:, 1), out_profile(:, 6), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 7), 'b');
hold off;

figure(7);
plot(IMU_true(:, 1), out_profile(:, 1), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 2), 'b');
hold off;

figure(8);
plot(IMU_true(:, 1), out_profile(:, 2), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 3), 'b');
hold off;

figure(9);
plot(IMU_true(:, 1), out_profile(:, 3), 'r');
hold on;
plot(IMU_true(:, 1), in_profile(2:end, 4), 'b');
hold off;



