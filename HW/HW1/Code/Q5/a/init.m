%% init
load('DataSet_02.mat');
accel = [IMU_true(:, 1), IMU_true(:, 2:4)];
gyro = [IMU_true(:, 1), IMU_true(:, 5:7)];