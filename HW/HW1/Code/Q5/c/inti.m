load('DataSet_02.mat')
Fb = timeseries(IMU_true(:, 2:4), IMU_true(:, 1));
wib = timeseries(IMU_true(:, 5:7), IMU_true(:, 1));