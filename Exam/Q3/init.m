%% init
load('DataSet_02.mat');
load('Q3.mat')
accel = [IMU_true(:, 1), IMU_true(:, 2:4)];
gyro = [IMU_true(:, 1), IMU_true(:, 5:7)];
dt = 0.01;
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
u_old = 0;
x_old = 0;
delta_x_num = zeros(9, 1);
%% run simulation
for i = 1:length(IMU_true)
    u_num = IMU_true(i, 2:7)';
    A_num = A_func(x_num, u_num);
    B_num = B_func(x_num, u_num);
    delta_u_num = u_num - u_old;
    delta_x_dot_num = A_num*delta_x_num + B_num*delta_u_num;
    delta_x_num = delta_x_num + delta_x_dot_num*dt;
    x_num = x_num + delta_x_num;
    u_old = u_num;
    out_profile(i, :) = x_num';
end